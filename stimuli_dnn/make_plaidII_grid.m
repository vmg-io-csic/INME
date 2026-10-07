function [configs, geo] = make_plaidII_grid(ioc_targets_deg, geometries, families)
% MAKE_PLAIDII_GRID  Config table crossing IOC directions with plaid geometries.
%
%   [configs, geo] = make_plaidII_grid(ioc_targets_deg, geometries, families)
%
%   Stacks one make_plaidII_configs sweep per geometry, giving a two-factor
%   design:
%
%     * DIRECTION (the rows within one geometry) — a rigid rotation of the
%       same plaid, so anything that changes across it is a property of the
%       model, not of the stimulus.
%     * GEOMETRY (the blocks) — a different plaid shape at every direction.
%       The speed ratio is the interesting axis here: it sets how far the
%       vector-average prediction falls from the IOC prediction, so a model
%       that stays on IOC as the ratio grows is passing a progressively
%       harder test, not just repeating the same one.
%
%   The two factors are crossed, so the direction effect and the geometry
%   effect can be read separately instead of being confounded.
%
% Inputs
%   ioc_targets_deg (1xK double) IOC directions to test [deg].
%                    Default -30:10:30.
%   geometries      (Gx2 double) one row per plaid shape:
%                    [speed_ratio, comp_sep_deg].  speed_ratio is
%                    speed2/speed1 (> 1); comp_sep_deg is the angular
%                    separation between the two component normals [deg].
%                    Default [1.5 20; 2 20; 3 20; 2 15; 2 30].
%   families        (char) 'A' | 'B' | 'both'.  Default 'both'.
%
% Outputs
%   configs (Cx5 double) [speed1, speed2, theta1_deg, theta2_deg, ioc_deg],
%                    C = K * G * numel(families).
%   geo     (struct) per-row diagnostics, same fields as make_plaidII_configs
%                    plus, all Cx1:
%             .speed_ratio   speed ratio of the row's geometry
%             .comp_sep_deg  component separation of the row's geometry
%             .geom_idx      which row of `geometries` this came from
%             .delta1_deg, .delta2_deg  component offsets from the pattern
%                            direction [deg]
%           and, 1xG summaries:
%             .geom_table    the `geometries` input, echoed
%             .geom_va_gap   |VA - IOC| for each geometry [deg] — the
%                            headroom the Type-II test has in that block
%
% See also: make_plaidII_configs, plaidII_predictions, run_plaidII_grid_2D

    if nargin < 1 || isempty(ioc_targets_deg), ioc_targets_deg = -30:10:30; end
    if nargin < 2 || isempty(geometries)
        geometries = [1.5 20; 2.0 20; 3.0 20; 2.0 15; 2.0 30];
    end
    if nargin < 3 || isempty(families), families = 'both'; end

    validateattributes(geometries, {'numeric'}, {'2d','ncols',2,'real'});

    n_geom = size(geometries, 1);
    n_fam  = 1 + strcmpi(families, 'both');
    n_g    = numel(ioc_targets_deg) * n_fam;   % rows contributed per geometry
    C      = n_geom * n_g;

    configs = zeros(C, 5);
    blank   = zeros(C, 1);
    geo = struct('ioc_deg', blank, 'va_deg', blank, 'va_ioc_diff', blank, ...
                 'family', blank, 'type', blank, 'speed_ratio', blank, ...
                 'comp_sep_deg', blank, 'geom_idx', blank, ...
                 'delta1_deg', blank, 'delta2_deg', blank);
    geom_va_gap = zeros(1, n_geom);

    for g = 1:n_geom
        % make_plaidII_configs validates the geometry and errors if the
        % requested ratio/separation would not give a Type-II plaid.
        [cfg_g, geo_g] = make_plaidII_configs(ioc_targets_deg, ...
            geometries(g,1), geometries(g,2), families);

        rows = (g-1)*n_g + (1:n_g);
        configs(rows, :) = cfg_g;

        geo.ioc_deg(rows)      = geo_g.ioc_deg;
        geo.va_deg(rows)       = geo_g.va_deg;
        geo.va_ioc_diff(rows)  = geo_g.va_ioc_diff;
        geo.family(rows)       = geo_g.family;
        geo.type(rows)         = geo_g.type;
        geo.speed_ratio(rows)  = geometries(g,1);
        geo.comp_sep_deg(rows) = geometries(g,2);
        geo.geom_idx(rows)     = g;
        geo.delta1_deg(rows)   = geo_g.delta_deg(1);
        geo.delta2_deg(rows)   = geo_g.delta_deg(2);

        geom_va_gap(g) = abs(geo_g.va_ioc_diff(1));
    end

    geo.geom_table  = geometries;
    geo.geom_va_gap = geom_va_gap;
end
