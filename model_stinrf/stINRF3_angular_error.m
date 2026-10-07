function err = stINRF3_angular_error(predicted_deg, target_deg)
% STINRF3_ANGULAR_ERROR  Signed angular difference between two directions.
%
% Inputs
%   predicted_deg (double, scalar or array) predicted direction(s) [deg]
%   target_deg (double, same size) expected direction(s) [deg]
%
% Outputs
%   err (double, same size) signed angular error in [-180, +180] deg

    err = mod(predicted_deg - target_deg + 180, 360) - 180;
end
