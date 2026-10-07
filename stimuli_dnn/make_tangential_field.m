function proj_raw = make_tangential_field(O_H, O_V, cx_wheel, cy_wheel)
% MAKE_TANGENTIAL_FIELD  per-pixel tangential projection
%   of (O_H,O_V) around (cx_wheel,cy_wheel)
%
%     proj(x) = OH(x)*e_tan_x(x) + OV(x)*e_tan_y(x)
%     e_tan = (1/r)*(-dy, dx),  dy = row-cy_wheel, dx = col-cx_wheel
%
%   Positive = clockwise (CW), negative = counter-clockwise (CCW)
%
%   proj_raw = make_tangential_field(O_H, O_V, cx_wheel, cy_wheel)
%
% Inputs
%   O_H, O_V (double, rows x cols) Tb-integrated motion energy maps
%   cx_wheel, cy_wheel (double) wheel centre [px]
%
% Outputs
%   proj_raw (double, rows x cols) unsmoothed tangential projection field

    [rows, cols] = size(O_H);
    [Xg, Yg] = meshgrid(1:cols, 1:rows);
    dx = Xg - cx_wheel;
    dy = Yg - cy_wheel;
    r_safe = max(sqrt(dx.^2 + dy.^2), 1e-6);

    Tx = -dy ./ r_safe;
    Ty =  dx ./ r_safe;
    proj_raw = O_H .* Tx + O_V .* Ty;
end
