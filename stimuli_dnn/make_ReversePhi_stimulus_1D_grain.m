function mI = make_ReversePhi_stimulus_1D_grain(N, n_frames, framerate, dov, ...
                sigma_t_ms, phase_mov, speed_dps, meanL, Type, Contrast, grain_px, seed)
% MAKE_REVERSEPHI_STIMULUS_1D_GRAIN  Reverse-phi with coarse-grained noise.
%
%   Same as make_ReversePhi_stimulus_1D.m but the binary noise is built in
%   blocks of grain_px pixels instead of one pixel per sample.
%
%   WHY. Reverse-phi needs correspondence between consecutive frames, and
%   white noise of 1 px grain is uncorrelated beyond 1 px: at the 20
%   px/frame displacement these networks need, there is nothing to match.
%   The networks also downsample internally by x8, which erases 1 px
%   structure outright. The rule is grain >= 2 * displacement, i.e. the
%   quarter-period rule applied to a noise field whose dominant period is
%   2*grain.
%
%   The frame-to-frame contrast inversion is untouched: it still flips
%   every time the integer pixel shift changes, which at these speeds is
%   every frame.
%
% Inputs
%   grain_px (int) width of each noise block [px]; 1 reproduces the
%                  original generator exactly
%   (the rest as make_ReversePhi_stimulus_1D.m)
%
% Outputs
%   mI (double, n_frames x N)

    if nargin < 11 || isempty(grain_px), grain_px = 1; end
    if nargin < 12, seed = []; end

    dx    = dov / N;                     % pixel size [deg]
    t_vec = 0 : n_frames-1;
    sigma_fr = sigma_t_ms / 1000 * framerate;
    G_env = exp(-(t_vec - n_frames/2).^2 / (2 * sigma_fr^2));

    mI = zeros(n_frames, N);

    % --- coarse-grained binary noise ---------------------------------
    % Without a seed the noise sample changes on every call and cannot be
    % controlled from outside. Passing `seed` reproduces one particular noise
    % realisation, which is what is needed to compare conditions without the
    % noise sample confounding them.
    if isempty(seed), rng('shuffle'); else, rng(seed); end
    n_blocks   = ceil(N / grain_px);
    block_vals = sign(randn(1, n_blocks));
    block_vals(block_vals == 0) = 1;
    base_noise = repelem(block_vals, grain_px);
    base_noise = base_noise(1:N);

    dt = 1 / framerate;
    current_contrast_sign = 1;
    last_pixel_shift = 0;

    for tt = 1:n_frames
        spatial_shift_deg = speed_dps * (tt * dt);
        pixel_shift = round(spatial_shift_deg / dx);

        if Type == -1
            if pixel_shift ~= last_pixel_shift
                current_contrast_sign = -current_contrast_sign;
                last_pixel_shift = pixel_shift;
            end
        else
            current_contrast_sign = 1;
        end

        shifted_noise = circshift(base_noise, pixel_shift);
        mI(tt, :) = meanL * (1 + Contrast * shifted_noise * ...
                             current_contrast_sign * G_env(tt));
    end
end
