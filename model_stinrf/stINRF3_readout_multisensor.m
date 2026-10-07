function [direction_deg, magnitude, info] = stINRF3_readout_multisensor(O_H, O_V, half_region)
% STINRF3_READOUT_MULTISENSOR  Extract direction from a square ROI mean.
%
% Inputs
%   O_H (double, rows x cols) horizontal channel (Tb-integrated)
%   O_V (double, rows x cols) vertical channel (Tb-integrated)
%   half_region (int, scalar) half-side of the square ROI [px]
%                ROI spans [centre-half_region+1 : centre+half_region].
%
% Outputs
%   direction_deg (double, scalar) motion direction [deg]
%   magnitude (double, scalar) sqrt(mean_OH^2 + mean_OV^2)
%   info (struct) diagnostic fields:
%                    .mean_OH, .mean_OV — ROI-averaged responses
%                    .all_directions    — per-pixel directions in the ROI
%                    .circ_std_deg      — circular std of per-pixel directions
%                    .R                 — mean resultant length
%                    .roi_bounds        — [r1, r2, c1, c2] pixel indices

    [rows, cols] = size(O_H);
    cx = round(cols / 2);
    cy = round(rows / 2);

    r1 = max(1, cy - half_region + 1);
    r2 = min(rows, cy + half_region);
    c1 = max(1, cx - half_region + 1);
    c2 = min(cols, cx + half_region);

    roi_H = O_H(r1:r2, c1:c2);
    roi_V = O_V(r1:r2, c1:c2);

    mean_OH = mean(roi_H(:));
    mean_OV = mean(roi_V(:));

    direction_deg = atan2d(mean_OV, mean_OH);
    magnitude     = sqrt(mean_OH^2 + mean_OV^2);

    all_dirs = atan2d(roi_V, roi_H);
    phi_rad  = deg2rad(all_dirs(:));
    R        = abs(mean(exp(1i * phi_rad)));

    info.mean_OH        = mean_OH;
    info.mean_OV        = mean_OV;
    info.all_directions = all_dirs;
    info.circ_std_deg   = stINRF3_circ_std(all_dirs(:));
    info.R              = R;
    info.roi_bounds     = [r1, r2, c1, c2];
end
