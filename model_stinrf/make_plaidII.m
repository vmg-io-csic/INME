function TCB = make_plaidII(N, n_frames, f_cpx, speed1, speed2, ...
                            theta1_deg, theta2_deg, phase1_rad, phase2_rad, ...
                            contrast, meanL)
% MAKE_PLAIDII  Generate a Type-II additive plaid stimulus.
%
% Inputs
%   N (int) image size [px] (square: N x N)
%   n_frames (int) number of temporal frames
%   f_cpx (double) spatial frequency [cycles/pixel] — same for
%                          both components.
%   speed1 (double) temporal frequency of component 1 [cycles/frame]
%   speed2 (double) temporal frequency of component 2 [cycles/frame]
%                          In the run script, these are computed from the
%                          config table as configs(c,1)*tf_scale and
%                          configs(c,2)*tf_scale, where tf_scale converts
%                          relative speed units to cy/frame.
%   theta1_deg (double) direction of component 1 [degrees]
%   theta2_deg (double) direction of component 2 [degrees]
%   phase1_rad (double) initial phase of component 1 [radians]
%   phase2_rad (double) initial phase of component 2 [radians]
%   contrast (double) michelson contrast per component, typically 1
%   meanL (double) mean luminance [a.u.], typically 100
%
% Outputs
%   TCB (single, N x N x n_frames) stimulus tensor

    th1 = deg2rad(theta1_deg);
    th2 = deg2rad(theta2_deg);
    coords = linspace(-N/2, N/2, N);
    [X, Y] = meshgrid(coords, coords);
    Xth1   = -X*cos(th1) + Y*sin(th1);
    Xth2   = -X*cos(th2) + Y*sin(th2);

    TCB = zeros(N, N, n_frames, 'single');
    for t = 1:n_frames
        g1 = contrast * sin(2*pi * f_cpx * Xth1 + phase1_rad + 2*pi*speed1*(t-1));
        g2 = contrast * sin(2*pi * f_cpx * Xth2 + phase2_rad + 2*pi*speed2*(t-1));
        TCB(:,:,t) = single(meanL * (1 + g1 + g2));
    end
end
