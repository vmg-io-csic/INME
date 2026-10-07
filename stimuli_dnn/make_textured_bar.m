function [TCB, info] = make_textured_bar(N, n_frames, v_px, direction_deg, ...
        bar_width_px, contrast, meanL, period_lo_px, period_hi_px, ...
        seed_bar, seed_bg, bg_v_px, edge_soft_px)
% MAKE_TEXTURED_BAR  A textured bar translating over a differently-textured
% background.
%
%   [TCB, info] = make_textured_bar(N, n_frames, v_px, direction_deg,
%           bar_width_px, contrast, meanL, period_lo_px, period_hi_px,
%           seed_bar, seed_bg, bg_v_px, edge_soft_px)
%
%   WHY THIS ONE
%   ------------
%   The textured grating (make_textured_grating) removes the aperture
%   problem but keeps a single global motion: every pixel moves the same
%   way. Natural training footage does not look like that. In Sintel almost
%   every frame contains a figure moving over a differently-textured
%   background, with a real occlusion boundary between them, and the
%   networks are trained to find exactly that boundary.
%
%   So this stimulus adds the one thing the co-moving texture leaves out:
%   two regions with different texture and different velocity, separated by
%   an edge that moves with the figure. If the networks recover the
%   direction here and not on the bare grating, the missing ingredient was
%   trackable 2D structure; if they recover it here and not on the textured
%   grating either, the missing ingredient was the figure/ground boundary.
%
%   GEOMETRY
%   --------
%   The bar's long axis is perpendicular to the motion direction and it
%   travels along its own normal, so the geometry matches the grating's
%   exactly and the same seven directions can be run side by side. Position
%   along the motion axis is p = x*cos(theta) - y*sin(theta), which is the
%   grating's projected coordinate negated, so p grows in the direction the
%   grating drifts.
%
%   Both textures are band-limited and exactly N-periodic, and both are
%   translated by a Fourier phase ramp rather than by interpolation, for
%   the reasons set out in make_textured_grating. The bar edge is softened
%   over a few pixels: a hard edge on a discrete grid is aliased, and the
%   point of this stimulus is the boundary, so the boundary should not be
%   the one part of the stimulus carrying a sampling artefact.
%
%   The bar wraps with period N. Its texture, the background texture and
%   the wrap are all N-periodic, so a frame is a rigid transformation of
%   the previous one everywhere, edges included.
%
%   One consequence of that wrap: the projection p spans +-N/2 at 0 and 90
%   degrees but +-N/2*(cos45+sin45) = +-0.71*N at 45, so obliquely the
%   canvas covers more than one wrap period and two bars are visible
%   instead of one. This is the same geometry the grating already has -- a
%   0.5 cpd grating shows 5.6 cycles across the canvas at 0 degrees and 8
%   at 45 -- and it is harmless here because the repeat period N is far
%   above twice the displacement, so no correspondence is ambiguous. It
%   does mean the stimulus is a periodic train of textured bands rather
%   than a single isolated object, which is if anything preferable: the
%   readout region always contains both textures and the figure/ground
%   boundary, at every direction and every frame.
%
% Inputs
%   N (int) canvas size in px, square N x N
%   n_frames (int) number of frames
%   v_px (double) bar speed [px/frame]
%   direction_deg (double) motion direction, 0 = rightward, 90 = upward
%   bar_width_px (double) bar width along the motion axis
%   contrast (double) RMS contrast of both textures
%   meanL (double) mean luminance
%   period_lo_px, period_hi_px (double) texture band, px per cycle
%   seed_bar, seed_bg (int) RNG seeds; must differ, that is the point
%   bg_v_px (double) background speed [px/frame]; 0 = static background
%   edge_soft_px (double) width of the bar's edge ramp, default 3
%
% Outputs
%   TCB (single, N x N x n_frames) luminance tensor
%   info (struct) .v_px_per_frame .vx_col .vy_row .bar_width_px
%                 .bg_v_px .peak_excursion .clipped .duty_cycle

    if nargin < 8  || isempty(period_lo_px),  period_lo_px = 40;  end
    if nargin < 9  || isempty(period_hi_px),  period_hi_px = 480; end
    if nargin < 10 || isempty(seed_bar),      seed_bar = 2;       end
    if nargin < 11 || isempty(seed_bg),       seed_bg  = 1;       end
    if nargin < 12 || isempty(bg_v_px),       bg_v_px  = 0;       end
    if nargin < 13 || isempty(edge_soft_px),  edge_soft_px = 3;   end
    if seed_bar == seed_bg
        error('make_textured_bar: seed_bar and seed_bg must differ');
    end

    theta  = deg2rad(direction_deg);
    coords = linspace(-N/2, N/2, N);
    [X, Y] = meshgrid(coords, coords);
    P = X*cos(theta) - Y*sin(theta);      % position along the motion axis

    vx_col =  v_px * cos(theta);
    vy_row = -v_px * sin(theta);

    T_bar = band_noise(N, period_lo_px, period_hi_px, seed_bar);
    T_bg  = band_noise(N, period_lo_px, period_hi_px, seed_bg);
    F_bar = fft2(T_bar);  F_bg = fft2(T_bg);

    k = 0:N-1;  k(k > N/2) = k(k > N/2) - N;
    [KX, KY] = meshgrid(k, k);
    FX = KX / N;  FY = KY / N;

    TCB = zeros(N, N, n_frames, 'single');
    duty = zeros(1, n_frames);
    for t = 1:n_frames
        dt = t - 1;

        ramp_bar = exp(-2i*pi * (FY*(vy_row*dt) + FX*(vx_col*dt)));
        Tb = real(ifft2(F_bar .* ramp_bar));
        if bg_v_px == 0
            Tg = T_bg;
        else
            r = bg_v_px / max(v_px, eps);
            ramp_bg = exp(-2i*pi * (FY*(vy_row*r*dt) + FX*(vx_col*r*dt)));
            Tg = real(ifft2(F_bg .* ramp_bg));
        end

        % bar mask: distance to the bar centre along the motion axis,
        % wrapped with period N so the bar re-enters seamlessly
        d = mod(P - v_px*dt + N/2, N) - N/2;
        M = smooth_band(abs(d), bar_width_px/2, edge_soft_px);
        duty(t) = mean(M(:));

        TCB(:,:,t) = single(meanL * (1 + contrast * (M.*Tb + (1-M).*Tg)));
    end

    lo = min(TCB(:));  hi = max(TCB(:));
    info = struct('v_px_per_frame', v_px, 'vx_col', vx_col, 'vy_row', vy_row, ...
                  'bar_width_px', bar_width_px, 'bg_v_px', bg_v_px, ...
                  'contrast', contrast, ...
                  'peak_excursion', max(abs([lo hi] - meanL))/meanL, ...
                  'clipped', lo < 0 || hi > 2*meanL, ...
                  'duty_cycle', mean(duty), ...
                  'period_lo_px', period_lo_px, 'period_hi_px', period_hi_px, ...
                  'seed_bar', seed_bar, 'seed_bg', seed_bg);
end


function T = band_noise(N, period_lo_px, period_hi_px, seed)
% Zero-mean, unit-RMS, exactly N-periodic noise confined to a radial band.
% No clipping: a pointwise clip would put energy back outside the band.
    rs = RandStream('mt19937ar', 'Seed', seed);
    F  = fft2(randn(rs, N, N));
    k  = 0:N-1;  k(k > N/2) = k(k > N/2) - N;
    [KX, KY] = meshgrid(k, k);
    R  = sqrt(KX.^2 + KY.^2) / N;
    F(~((R >= 1/period_hi_px) & (R <= 1/period_lo_px))) = 0;
    T = real(ifft2(F));
    T = T - mean(T(:));
    T = T / std(T(:));
end


function M = smooth_band(dist, half_width, soft)
% 1 inside the band, 0 outside, with a smoothstep ramp of width `soft`
% across the edge so the boundary is not aliased by the pixel grid.
    u = (half_width + soft/2 - dist) / max(soft, eps);
    u = max(min(u, 1), 0);
    M = u.*u.*(3 - 2*u);
end
