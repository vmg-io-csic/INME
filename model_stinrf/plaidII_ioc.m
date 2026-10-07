function [ioc_deg, va_deg, ioc_speed, is_typeII] = plaidII_ioc(speed1, speed2, theta1_deg, theta2_deg)
% PLAIDII_IOC  Intersection-of-Constraints and Vector-Average predictions.
%
% Inputs
%   speed1 (double) relative normal speed of component 1
%   speed2 (double) relative normal speed of component 2
%   theta1_deg (double) direction of component 1 [deg]
%   theta2_deg (double) direction of component 2 [deg]
%
% Outputs
%   ioc_deg (double) IOC pattern direction [deg]
%   va_deg (double) VA pattern direction [deg]
%   ioc_speed (double) IOC pattern speed, same units as speed1/speed2
%   is_typeII (logical) true when the IOC direction falls OUTSIDE the
%                     angular wedge spanned by the two components, which is
%                     the defining property of a Type-II plaid

    A = [cosd(theta1_deg), sind(theta1_deg); ...
         cosd(theta2_deg), sind(theta2_deg)];

    if abs(det(A)) < 1e-12
        error('plaidII_ioc:parallel', ...
              'Components are parallel (theta1 = %g, theta2 = %g): IOC undefined.', ...
              theta1_deg, theta2_deg);
    end

    v = A \ [speed1; speed2];
    ioc_deg   = atan2d(v(2), v(1));
    ioc_speed = hypot(v(1), v(2));

    va = [speed1*cosd(theta1_deg) + speed2*cosd(theta2_deg); ...
          speed1*sind(theta1_deg) + speed2*sind(theta2_deg)] / 2;
    va_deg = atan2d(va(2), va(1));

    % Type-II test: is the IOC direction inside the wedge theta1 -> theta2 ?
    span   = wrap180(theta2_deg - theta1_deg);
    offset = wrap180(ioc_deg    - theta1_deg);
    if span >= 0
        inside = (offset >= 0) && (offset <= span);
    else
        inside = (offset <= 0) && (offset >= span);
    end
    is_typeII = ~inside;
end

function a = wrap180(a)
    a = mod(a + 180, 360) - 180;
end
