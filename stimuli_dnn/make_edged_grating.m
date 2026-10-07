function [TCB, info] = make_edged_grating(N, n_frames, f_cpx, speed, ...
        direction_deg, phase0_rad, contrast, meanL, bar_width_px, ...
        edge_soft_px)
% MAKE_EDGED_GRATING  A drifting grating confined to a moving band, on a
% uniform background.
%
%   [TCB, info] = make_edged_grating(N, n_frames, f_cpx, speed,
%           direction_deg, phase0_rad, contrast, meanL, bar_width_px,
%           edge_soft_px)
%
%   WHY THIS ONE
%   ------------
%   It fills the one empty cell of a 2x2. The other three exist already:
%
%                     no border                 with border
%     no texture      bare grating              THIS
%     texture         make_textured_grating     make_textured_bar
%
%   The textured bar changes two things at once relative to the bare
%   grating: it adds trackable 2D texture and it adds a figure/ground
%   boundary. Dual and ME-Attention recover the direction from it but not
%   from the textured grating, so for them the boundary looks like the
%   ingredient that matters -- but that cannot be concluded while no
%   stimulus has a boundary and no texture. This one does: the carrier is
%   the same sinusoid as ever, and the only thing added is the edge.
%
%   WHY THE WINDOW MOVES
%   --------------------
%   With a stationary aperture, the line terminators at the edge slide
%   along it as the grating drifts and carry no information about the
%   pattern's true direction -- that is the classical aperture-problem
%   demonstration, and such a stimulus would be no more informative than
%   the bare grating. Here the window translates with the grating, so the
%   terminators are genuine 2D features moving at the true velocity, and a
%   system that tracks them recovers the direction.
%
%   GEOMETRY
%   --------
%   The band's long axis is perpendicular to the motion direction and it
%   travels along its own normal, matching make_textured_bar exactly so
%   the two can be compared width for width. Position along the motion
%   axis is p = x*cos(theta) - y*sin(theta); the band wraps with period N,
%   so obliquely the canvas shows more than one band, as it does there.
%
%   Outside the band the field is flat at meanL: no texture is added
%   anywhere, which is the whole point.
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
%   bar_width_px (double) band width along the motion axis
%   edge_soft_px (double) width of the edge ramp, default 3
%
% Outputs
%   TCB (single, N x N x n_frames) luminance tensor
%   info (struct) .v_px_per_frame .vx_col .vy_row .bar_width_px
%                 .duty_cycle .max_excursion .clipped

    if nargin < 10 || isempty(edge_soft_px), edge_soft_px = 3; end

    theta  = deg2rad(direction_deg);
    coords = linspace(-N/2, N/2, N);
    [X, Y] = meshgrid(coords, coords);
    P = X*cos(theta) - Y*sin(theta);       % along the motion axis
    X_th = -P;                             % make_grating's projected coord

    v_px   = speed / f_cpx;
    vx_col =  v_px * cos(theta);
    vy_row = -v_px * sin(theta);

    TCB = zeros(N, N, n_frames, 'single');
    duty = zeros(1, n_frames);
    for t = 1:n_frames
        dt = t - 1;
        G = contrast * sin(2*pi * f_cpx * X_th + phase0_rad + 2*pi*speed*dt);
        d = mod(P - v_px*dt + N/2, N) - N/2;
        M = smooth_band(abs(d), bar_width_px/2, edge_soft_px);
        duty(t) = mean(M(:));
        TCB(:,:,t) = single(meanL * (1 + M .* G));
    end

    lo = min(TCB(:));  hi = max(TCB(:));
    info = struct('v_px_per_frame', v_px, 'vx_col', vx_col, 'vy_row', vy_row, ...
                  'bar_width_px', bar_width_px, 'contrast', contrast, ...
                  'duty_cycle', mean(duty), ...
                  'max_excursion', max(abs([lo hi] - meanL))/meanL, ...
                  'clipped', lo < 0 || hi > 2*meanL);
end


function M = smooth_band(dist, half_width, soft)
% 1 inside the band, 0 outside, smoothstep across the edge so the boundary
% is not aliased by the pixel grid.
    u = (half_width + soft/2 - dist) / max(soft, eps);
    u = max(min(u, 1), 0);
    M = u.*u.*(3 - 2*u);
end
