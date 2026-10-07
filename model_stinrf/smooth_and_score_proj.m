function [proj_tan, proj_field] = smooth_and_score_proj(proj_raw, ...
                                cx_wheel, cy_wheel, cx_canvas, cy_canvas, ...
                                r_inner, r_outer, sigma_profile_name)
% SMOOTH_AND_SCORE_PROJ
%
%   [proj_tan, proj_field] = smooth_and_score_proj(proj_raw, ...
%       cx_wheel, cy_wheel, cx_canvas, cy_canvas, r_inner, r_outer, ...
%       sigma_profile_name)
%
%
% Inputs
%   proj_raw (double, rows x cols) unsmoothed tangential projection field
%             (from make_tangential_field.m)
%   cx_wheel, cy_wheel (double) wheel centre [px]
%   cx_canvas, cy_canvas (double) canvas centre / fovea [px]
%   r_inner, r_outer (double) annulus radii around the wheel centre [px]
%   sigma_profile_name (char, optional) 'A','C','D','F','Unif20','Unif40'
%             (default 'D'): exponential smoothing profile D
%
% Outputs
%   proj_tan (double) global scalar, positive = CW
%   proj_field (double, rows x cols) the smoothed per-pixel field

    if nargin < 8 || isempty(sigma_profile_name), sigma_profile_name = 'D'; end

    [rows, cols] = size(proj_raw);
    [Xg, Yg] = meshgrid(1:cols, 1:rows);

    dxc = Xg - cx_canvas;
    dyc = Yg - cy_canvas;
    rc  = sqrt(dxc.^2 + dyc.^2);
    sigma_map = local_sigma_profile(rc, sigma_profile_name);
    proj_field = variable_gauss(proj_raw, sigma_map);

    dxw = Xg - cx_wheel;
    dyw = Yg - cy_wheel;
    rw  = sqrt(dxw.^2 + dyw.^2);
    mask = (rw >= r_inner) & (rw < r_outer);
    if ~any(mask(:))
        proj_tan = NaN;
    else
        proj_tan = mean(proj_field(mask));
    end
end

function sigma = local_sigma_profile(r, name)
% LOCAL_SIGMA_PROFILE
    switch name
        case 'A'        % linear 78 -> 20
            r_max = max(r(:));
            sigma = 78.0 - 58.0 * (r / r_max);
        case 'C'        % linear 40 -> 5
            r_max = max(r(:));
            sigma = 40.0 - 35.0 * (r / r_max);
        case 'D'        % exponential 80 -> 5 (validated optimum)
            sigma = 5.0 + 75.0 * exp(-r / 30.0);
        case 'F'        % exponential 40 -> 5
            sigma = 5.0 + 35.0 * exp(-r / 40.0);
        case 'Unif20'
            sigma = 20.0 * ones(size(r));
        case 'Unif40'
            sigma = 40.0 * ones(size(r));
        otherwise
            error('smooth_and_score_proj: unknown sigma_profile_name ''%s''.', name);
    end
end
