function s = stINRF3_circ_std(angles_deg)
% STINRF3_CIRC_STD  Circular standard deviation of a set of angles.
%
% Inputs
%   angles_deg (double, vector) angles in degrees
%
% Outputs
%   s (double, scalar) circular standard deviation in degrees

    phi   = deg2rad(angles_deg(:));
    R     = abs(mean(exp(1i * phi)));
    s_raw = sqrt(-2 * log(R + eps));
    if isreal(s_raw)
        s = rad2deg(s_raw);
    else
        s = std(angles_deg(:));
    end
end
