function [configs, geo] = make_plaidII_configs(ioc_targets_deg, speed_ratio, comp_sep_deg, families)
% MAKE_PLAIDII_CONFIGS  Build a Type-II plaid config table by inverse design.
%
%   [configs, geo] = make_plaidII_configs(ioc_targets_deg, speed_ratio, comp_sep_deg, families)
%
%   Instead of hand-picking component directions and reading off the
%   resulting IOC direction, this solves the inverse problem: given the
%   IOC directions you want to test, it returns the plaids that produce
%   them exactly.
%
%   GEOMETRY.  Write the two component normals as offsets from the pattern
%   (IOC) direction: delta_i = (component normal) - (IOC direction).  The
%   IOC constraint V.n_i = v_i then reads
%
%       speed_i = |V| * f * cos(delta_i),
%
%   so the SPEED RATIO fixes the two offsets and the IOC direction only
%   fixes their common rotation:
%
%       r = speed2/speed1 = cos(delta2)/cos(delta1),   delta2 = delta1 - sep
%       =>  tan(delta1) = (r - cos sep) / sin sep
%
%   Consequences that make this the right way to sweep the phenomenon:
%     * Rotating both components rigidly rotates the IOC direction by the
%       same amount and leaves the speed ratio, the component separation
%       and the VA-IOC discrepancy untouched.  Every row of the sweep is
%       therefore the SAME Type-II plaid at a different orientation, so any
%       variation in the model output is a property of the model, not a
%       side effect of having changed the stimulus geometry too.
%     * delta1 and delta2 come out with the same sign (both components on
%       one side of the pattern motion), which is precisely the definition
%       of a Type-II plaid, and VA lands far from IOC.
%
%   The default speed_ratio = 2 and comp_sep_deg = 20 reproduce the
%   original three hand-made configs: they give delta = [72.122, 52.122],
%   i.e. exactly the geometry of [1 2 -45 -25 27.12] and [1 2 -50 -30 22.12].
%
%   FAMILIES.  For a given IOC direction the mirror-image plaid (both
%   components reflected across the pattern direction, delta -> -delta) has
%   the same IOC and the mirrored VA.  Running both families separates a
%   genuine dependence on pattern direction from a bias tied to the
%   components' own directions — the original table mixed the two families
%   (config 2 was the mirror of configs 1 and 3), which confounds them.
%
% Inputs
%   ioc_targets_deg (1xK double) IOC directions to test [deg], readout
%                    convention.  Default -30:5:30.
%   speed_ratio     (double) speed2/speed1, > 1.  Sets how far VA falls
%                    from IOC; 1 would give a Type-I plaid.  Default 2.
%   comp_sep_deg    (double) angular separation between the two component
%                    normals [deg], in (0, 180).  Default 20.
%   families        (char) 'A' | 'B' | 'both'.  'A' puts both components
%                    clockwise of the pattern direction, 'B' is the mirror
%                    image, 'both' stacks A then B.  Default 'both'.
%
% Outputs
%   configs (Cx5 double) one row per plaid, in the format the run scripts
%                    expect: [speed1, speed2, theta1_deg, theta2_deg, ioc_deg].
%                    speed1 is always 1 (relative units; the run scripts
%                    multiply by tf_scale).
%   geo     (struct) design diagnostics, one entry per row unless noted:
%             .ioc_deg      (Cx1) IOC direction, re-derived from the built
%                           row with plaidII_predictions (not copied)
%             .va_deg       (Cx1) vector-average direction [deg]
%             .va_ioc_diff  (Cx1) signed VA - IOC [deg]
%             .family       (Cx1) 1 for A, 2 for B
%             .type         (Cx1) plaid type as classified from the row
%             .delta_deg    (1x2) component offsets of family A [deg]
%             .speed_ratio, .comp_sep_deg  (scalars) echoed inputs
%
% See also: plaidII_predictions, make_plaidII, run_plaidII_2D

    if nargin < 1 || isempty(ioc_targets_deg), ioc_targets_deg = -30:5:30; end
    if nargin < 2 || isempty(speed_ratio),     speed_ratio     = 2;        end
    if nargin < 3 || isempty(comp_sep_deg),    comp_sep_deg    = 20;       end
    if nargin < 4 || isempty(families),        families        = 'both';   end

    validateattributes(speed_ratio,  {'numeric'}, {'scalar','real','positive'});
    validateattributes(comp_sep_deg, {'numeric'}, {'scalar','real','>',0,'<',180});

    if speed_ratio == 1
        error(['make_plaidII_configs: speed_ratio = 1 makes the two ' ...
               'components equal in speed, i.e. a Type-I plaid (VA = IOC). ' ...
               'Use a ratio > 1 (2 is the validated value).']);
    end

    % Component offsets from the pattern direction implied by the ratio
    delta1 = atand((speed_ratio - cosd(comp_sep_deg)) / sind(comp_sep_deg));
    delta2 = delta1 - comp_sep_deg;

    % A component with |delta| >= 90 would need a negative normal speed:
    % geometrically valid but no longer the intended stimulus.
    if abs(delta1) >= 90 || abs(delta2) >= 90
        error(['make_plaidII_configs: speed_ratio = %.3g with a %.3g' char(176) ...
               ' separation puts a component at %.2f' char(176) ' from the pattern ' ...
               'direction (must stay under 90). Lower the ratio or widen ' ...
               'the separation.'], speed_ratio, comp_sep_deg, max(abs([delta1 delta2])));
    end
    if sign(delta1) ~= sign(delta2)
        error(['make_plaidII_configs: the requested ratio/separation ' ...
               'straddles the pattern direction, which gives a Type-I ' ...
               'plaid (VA = IOC). Increase speed_ratio or reduce comp_sep_deg.']);
    end

    switch lower(families)
        case 'a',    signs = +1;
        case 'b',    signs = -1;
        case 'both', signs = [+1 -1];
        otherwise
            error('make_plaidII_configs: families must be ''A'', ''B'' or ''both''.');
    end

    ioc_targets_deg = ioc_targets_deg(:);
    K = numel(ioc_targets_deg);
    C = K * numel(signs);

    configs = zeros(C, 5);
    geo.ioc_deg     = zeros(C, 1);
    geo.va_deg      = zeros(C, 1);
    geo.va_ioc_diff = zeros(C, 1);
    geo.family      = zeros(C, 1);
    geo.type        = zeros(C, 1);

    row = 0;
    for k_fam = 1:numel(signs)
        sgn = signs(k_fam);
        for k = 1:K
            row = row + 1;
            ioc = ioc_targets_deg(k);

            % theta is the make_plaidII direction parameter; a component
            % drifting towards angle A is written with theta = -A, and
            % A = (IOC direction, math convention) + delta.
            theta1 = ioc - sgn*delta1;
            theta2 = ioc - sgn*delta2;

            configs(row, :) = [1, speed_ratio, theta1, theta2, ioc];

            % Re-derive from the row we just built, so the table is
            % verified by the forward solver rather than trusted.
            [ioc_chk, va_chk, in] = plaidII_predictions(1, speed_ratio, theta1, theta2);
            if abs(wrap180(ioc_chk - ioc)) > 1e-6
                error(['make_plaidII_configs: internal inconsistency at row %d ' ...
                       '(requested IOC %.4f, built %.4f).'], row, ioc, ioc_chk);
            end

            configs(row, 5)   = ioc_chk;
            geo.ioc_deg(row)  = ioc_chk;
            geo.va_deg(row)   = va_chk;
            geo.va_ioc_diff(row) = in.va_ioc_diff;
            geo.family(row)   = 1 + (sgn < 0);   % 1 = A, 2 = B (mirror)
            geo.type(row)     = in.type;
        end
    end

    if any(geo.type ~= 2)
        error('make_plaidII_configs: %d of %d rows are not Type-II plaids.', ...
              sum(geo.type ~= 2), C);
    end

    geo.delta_deg   = [delta1 delta2];
    geo.speed_ratio = speed_ratio;
    geo.comp_sep_deg = comp_sep_deg;
end

function a = wrap180(a)
    a = mod(a + 180, 360) - 180;
end
