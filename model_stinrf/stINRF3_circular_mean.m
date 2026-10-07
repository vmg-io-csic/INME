function [mean_dir, std_dir] = stINRF3_circular_mean(directions_deg)
% STINRF3_CIRCULAR_MEAN  Circular mean and standard deviation of directions.
%
% Inputs
%   directions_deg (double, vector) direction angles in degrees
%
% Outputs
%   mean_dir (double, scalar) circular mean direction in degrees
%   std_dir (double, scalar) circular standard deviation in degrees
%                                (via stINRF3_circ_std).

    phi      = deg2rad(directions_deg(:));
    z_mean   = mean(exp(1i * phi));
    mean_dir = rad2deg(angle(z_mean));
    std_dir  = stINRF3_circ_std(directions_deg(:));
end
