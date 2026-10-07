function [TCB, info] = make_slit_grating(N, n_frames, f_cpx, speed, ...
        direction_deg, phase0_rad, contrast, meanL, slit_width_px, ...
        edge_soft_px)
% MAKE_SLIT_GRATING  A drifting grating seen through a slit whose long axis
% runs along the direction of motion, so the bars are cut and their end
% points are visible.
%
%   [TCB, info] = make_slit_grating(N, n_frames, f_cpx, speed,
%           direction_deg, phase0_rad, contrast, meanL, slit_width_px,
%           edge_soft_px)
%
%   WHY THE SLIT RUNS ALONG THE MOTION
%   ----------------------------------
%   A grating's bars are lines of constant p, where p is the coordinate
%   along the motion. An aperture whose edges are also lines of constant p
%   -- a band perpendicular to the motion, which is what
%   make_edged_grating built -- cuts nothing: its edges are parallel to the
%   bars, every visible bar still runs the full length of the canvas, and
%   no end point exists anywhere. Such a stimulus adds a border but no
%   two-dimensional information at all, and the aperture problem survives
%   it untouched.
%
%   Turning the slit ninety degrees fixes that. Here the slit is bounded in
%   q, the coordinate perpendicular to the motion, so its two edges are
%   perpendicular to the bars. Every visible bar segment now runs from
%   q = -w/2 to q = +w/2 and ends at two terminators. As the grating
%   drifts, a terminator stays on its edge and slides along it at exactly
%   the pattern speed, so the terminators are genuine two-dimensional
%   features moving at the true velocity. This is the barber-pole geometry
%   with the pole aligned to the real motion, the case in which terminators
%   and pattern agree instead of competing.
%
%   The slit itself does not move: translating a strip along its own length
%   changes nothing in the image. The moving features are the bars and
%   their end points, which is the point.
%
%   Outside the slit the field is flat at meanL. No texture is added
%   anywhere, so this differs from the bare grating in exactly one respect:
%   the bars now have ends.
%
% Inputs
%   N (int) canvas size in px, square N x N
%   n_frames (int) number of frames
%   f_cpx (double) grating spatial frequency [cycles/px]
%   speed (double) grating temporal frequency [cycles/frame]
%   direction_deg (double) motion direction, 0 = rightward, 90 = upward
%   phase0_rad (double) initial spatial phase
%   contrast (double) grating Michelson contrast
%   meanL (double) mean luminance
%   slit_width_px (double) slit width across the motion axis
%   edge_soft_px (double) width of the edge ramp, default 3
%
% Outputs
%   TCB (single, N x N x n_frames) luminance tensor
%   info (struct) .v_px_per_frame .slit_width_px .duty_cycle
%                 .max_excursion .clipped

    if nargin < 10 || isempty(edge_soft_px), edge_soft_px = 3; end

    theta  = deg2rad(direction_deg);
    coords = linspace(-N/2, N/2, N);
    [X, Y] = meshgrid(coords, coords);
    P = X*cos(theta) - Y*sin(theta);       % along the motion
    Q = X*sin(theta) + Y*cos(theta);       % across it
    X_th = -P;                             % make_grating's projected coord

    v_px = speed / f_cpx;

    % Static slit, centred, bounded across the motion axis.
    M = smooth_band(abs(Q), slit_width_px/2, edge_soft_px);
    duty = mean(M(:));

    TCB = zeros(N, N, n_frames, 'single');
    for t = 1:n_frames
        G = contrast * sin(2*pi * f_cpx * X_th + phase0_rad + 2*pi*speed*(t-1));
        TCB(:,:,t) = single(meanL * (1 + M .* G));
    end

    lo = min(TCB(:));  hi = max(TCB(:));
    info = struct('v_px_per_frame', v_px, 'slit_width_px', slit_width_px, ...
                  'contrast', contrast, 'duty_cycle', duty, ...
                  'max_excursion', max(abs([lo hi] - meanL))/meanL, ...
                  'clipped', lo < 0 || hi > 2*meanL);
end


function M = smooth_band(dist, half_width, soft)
    u = (half_width + soft/2 - dist) / max(soft, eps);
    u = max(min(u, 1), 0);
    M = u.*u.*(3 - 2*u);
end
