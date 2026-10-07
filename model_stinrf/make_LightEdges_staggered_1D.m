function [stim_1D, lum_mean] = make_LightEdges_staggered_1D(N, n_frames, framerate, dov, ...
                                    phase0_rad, contrast, meanL, v_dps, period_deg, n_tiles)
% MAKE_LIGHTEDGES_STAGGERED_1D  Tessellated light edges with staggered resets.
%
%   Same light-edge stimulus as make_TimeReversal_stimulus_1D, except that
%   each spatial period is split into n_tiles tiles whose reset phases are
%   staggered by 1/n_tiles of a cycle. With n_tiles = 1 the tiles reset
%   simultaneously and the whole field steps from light to dark once per
%   cycle; with larger n_tiles the resets are spread across space and the
%   field-average luminance stays nearly constant.
%
% Inputs
%   N (int) number of spatial samples
%   n_frames (int) number of temporal frames
%   framerate (double) frames per second [fps]
%   dov (double) field of view [degrees]
%   phase0_rad (double) initial spatial phase [radians]
%   contrast (double) multiplier on the +/-1 profile
%   meanL (double) mean luminance [a.u.]
%   v_dps (double) edge speed [degrees/second]
%   period_deg (double) spatial period [degrees]
%   n_tiles (int) tiles per period; 1 reproduces the synchronised stimulus
%
% Outputs
%   stim_1D (double, n_frames x N) luminance
%   lum_mean (double, n_frames x 1) field-average of the +/-1 profile,
%       for diagnosing the size of the global luminance step

    x_deg   = linspace(-dov/2, dov/2, N);
    t_s     = (0:n_frames-1) / framerate;
    x_shift = phase0_rad * period_deg / (2*pi);

    u = mod(x_deg + x_shift, period_deg) / period_deg;   % position within period
    k = min(floor(u * n_tiles), n_tiles - 1);            % tile index
    u_tile = u * n_tiles - k;                            % position within tile

    c = zeros(n_frames, N);
    for tt = 1:n_frames
        edge = mod(v_dps * t_s(tt) / period_deg + k / n_tiles, 1);
        row  = -ones(1, N);
        row(u_tile <= edge) = 1;
        c(tt,:) = row;
    end

    lum_mean = mean(c, 2);
    stim_1D  = meanL * (1 + contrast * c);
end
