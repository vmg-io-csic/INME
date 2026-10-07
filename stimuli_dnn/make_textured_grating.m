function [TCB, info] = make_textured_grating(N, n_frames, f_cpx, speed, ...
        direction_deg, phase0_rad, c_grating, c_texture, meanL, ...
        period_lo_px, period_hi_px, seed)
% MAKE_TEXTURED_GRATING  Drifting grating plus a rigidly co-moving texture.
%
%   [TCB, info] = make_textured_grating(N, n_frames, f_cpx, speed,
%           direction_deg, phase0_rad, c_grating, c_texture, meanL,
%           period_lo_px, period_hi_px, seed)
%
%   WHY
%   ---
%   A drifting sinusoidal grating is a one-dimensional pattern: it has no
%   feature that can be tracked across frames, so only the component of
%   velocity normal to its bars is recoverable from the image sequence.
%   That is the aperture problem, and it is what the optical-flow networks
%   run into on the grating stimuli -- they return the correct normal speed
%   and an unconstrained tangential slide.
%
%   This generator adds two-dimensional structure that translates with
%   exactly the same velocity vector as the grating, so the composite has a
%   single, unambiguous 2D motion. If the networks recover the direction
%   here but not on the bare grating, the failure was the aperture problem
%   and nothing else.
%
%   HOW THE TEXTURE IS MOVED
%   ------------------------
%   By a phase ramp in the Fourier domain, not by interpolation. The shift
%   is 34.13 px/frame, not an integer, and a bilinear resample would blur
%   the texture a little more on every frame -- a slow loss of contrast
%   that the networks would see as something other than rigid motion. A
%   phase ramp translates an exactly band-limited, exactly N-periodic
%   texture by any real amount with no interpolation and no loss.
%
%   Because the texture is N-periodic, wrapping at the canvas edge is
%   seamless: the whole frame is a rigid translation of the previous one,
%   edges included. The grating is N-periodic in general it is not, but it
%   is generated analytically per frame, so it never wraps either.
%
%   CHOICE OF TEXTURE BAND
%   ----------------------
%   Two constraints, both from the 34.13 px/frame displacement:
%
%     - any texture component whose period is below twice the displacement
%       is itself ambiguous (a shift of d looks like a shift of d - lambda),
%       so periods must stay above ~68 px;
%     - structure finer than ~16 px does not survive the x8 internal
%       downsampling of these networks.
%
%   The default band [96, 480] px clears the first with margin and brackets
%   the grating's own period (170.7 px at 0.5 cpd, 85.33 px/deg).
%
% Inputs
%   N (int) canvas size in px, square N x N
%   n_frames (int) number of frames
%   f_cpx (double) grating spatial frequency [cycles/px]
%   speed (double) grating temporal frequency [cycles/frame]
%   direction_deg (double) motion direction, 0 = rightward, 90 = upward
%   phase0_rad (double) initial spatial phase of the grating
%   c_grating (double) grating Michelson contrast
%   c_texture (double) texture RMS contrast; 0 gives the bare grating
%   meanL (double) mean luminance
%   period_lo_px, period_hi_px (double) texture band, in px per cycle
%   seed (int) RNG seed for the texture
%
% Outputs
%   TCB (single, N x N x n_frames) luminance tensor
%   info (struct) .v_px_per_frame .vx_col .vy_row .texture_rms
%                 .max_excursion .clipped (true if the tensor hit 0 or 2*meanL)

    if nargin < 10 || isempty(period_lo_px), period_lo_px = 40;  end
    if nargin < 11 || isempty(period_hi_px), period_hi_px = 480; end
    if nargin < 12 || isempty(seed),         seed = 1;           end

    theta  = deg2rad(direction_deg);
    coords = linspace(-N/2, N/2, N);
    [X, Y] = meshgrid(coords, coords);
    X_th   = -X*cos(theta) + Y*sin(theta);

    % Velocity of the pattern, in array terms. Rows grow downwards while
    % the direction convention has 90 deg pointing up, hence the sign on
    % the row component.
    v_px   = speed / f_cpx;
    vx_col =  v_px * cos(theta);
    vy_row = -v_px * sin(theta);

    % ---- texture: band-limited, unit RMS, exactly N-periodic ----
    T0 = zeros(N, N);
    if c_texture > 0
        rs = RandStream('mt19937ar', 'Seed', seed);
        F  = fft2(randn(rs, N, N));
        k  = 0:N-1;  k(k > N/2) = k(k > N/2) - N;     % signed wavenumbers
        [KX, KY] = meshgrid(k, k);
        R  = sqrt(KX.^2 + KY.^2) / N;                 % cycles/px, radial
        keep = (R >= 1/period_hi_px) & (R <= 1/period_lo_px);
        F(~keep) = 0;
        T0 = real(ifft2(F));
        T0 = T0 - mean(T0(:));
        T0 = T0 / std(T0(:));
        % No clipping here. Clipping is a pointwise nonlinearity, and it
        % puts energy back at every spatial frequency -- including the fine
        % scales the band limit exists to remove. The band-limited field is
        % left exactly as it is, and the caller keeps the composite inside
        % the display range by choosing c_texture against info.peak_texture.
    end
    FT = fft2(T0);
    k  = 0:N-1;  k(k > N/2) = k(k > N/2) - N;
    [KX, KY] = meshgrid(k, k);
    FX = KX / N;  FY = KY / N;                        % cycles/px

    TCB = zeros(N, N, n_frames, 'single');
    for t = 1:n_frames
        dt = t - 1;
        phase_t = phase0_rad + 2*pi * speed * dt;
        G = c_grating * sin(2*pi * f_cpx * X_th + phase_t);
        if c_texture > 0
            % shift theorem: T(r - dr, c - dc) <-> F .* exp(-2i*pi*(fy*dr + fx*dc))
            ramp = exp(-2i*pi * (FY * (vy_row*dt) + FX * (vx_col*dt)));
            Tt = real(ifft2(FT .* ramp));
        else
            Tt = 0;
        end
        TCB(:,:,t) = single(meanL * (1 + G + c_texture * Tt));
    end

    lo = min(TCB(:));  hi = max(TCB(:));
    info = struct('v_px_per_frame', v_px, 'vx_col', vx_col, 'vy_row', vy_row, ...
                  'texture_rms', c_texture, 'peak_texture', max(abs(T0(:))), ...
                  'max_excursion', max(abs([lo hi] - meanL))/meanL, ...
                  'clipped', lo < 0 || hi > 2*meanL, ...
                  'period_lo_px', period_lo_px, 'period_hi_px', period_hi_px, ...
                  'seed', seed);
end
