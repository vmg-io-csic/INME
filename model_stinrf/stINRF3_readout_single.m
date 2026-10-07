function [direction_deg, magnitude] = stINRF3_readout_single(O_H, O_V)
% STINRF3_READOUT_SINGLE  Extract motion direction from the central pixel.
%
% Inputs
%   O_H (double, rows x cols) horizontal channel response (Tb-integrated)
%   O_V (double, rows x cols) vertical channel response (Tb-integrated)
%
% Outputs
%   direction_deg (double, scalar) motion direction in degrees
%                  Convention: 0° = rightward, 90° = upward.
%   magnitude (double, scalar) response magnitude = sqrt(O_H^2 + O_V^2)

    cx = round(size(O_H, 2) / 2);
    cy = round(size(O_H, 1) / 2);
    oh = O_H(cy, cx);
    ov = O_V(cy, cx);
    direction_deg = atan2d(ov, oh);
    magnitude     = sqrt(oh^2 + ov^2);
end
