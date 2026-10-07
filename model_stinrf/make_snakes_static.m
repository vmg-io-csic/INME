function TCB = make_snakes_static(N, n_frames, onset_f, offset_f, ...
                                   period_px, orientation_deg, L, meanL, ...
                                   phase_shift_px)
% MAKE_SNAKES_STATIC  Generate a static sawtooth pattern (appearance paradigm).
%
% Inputs
%   N (int) image size [px] (square: N x N)
%   n_frames (int) total number of frames
%   onset_f (int) first frame where the pattern is visible
%   offset_f (int) last frame where the pattern is visible
%   period_px (double) spatial period of the sawtooth [px]
%   orientation_deg (double) pattern orientation [degrees]
%   L (1x4 double) luminance levels per quarter-period
%   meanL (double) background luminance [a.u.]
%   phase_shift_px (double) spatial phase offset [px] (default: 0)
%
% Outputs
%   TCB (single, N x N x n_frames) stimulus tensor

    if nargin < 9, phase_shift_px = 0; end

    [X, Y] = meshgrid(0:N-1, 0:N-1);
    X = X - (N-1)/2;
    Y = Y - (N-1)/2;
    coord = -X*cosd(orientation_deg) + Y*sind(orientation_deg) + phase_shift_px;

    patch = build_sawtooth(coord, period_px, L, meanL, N);

    TCB = single(meanL * ones(N, N, n_frames));
    for t = onset_f : min(offset_f, n_frames)
        TCB(:,:,t) = patch;
    end
end

function patch = build_sawtooth(coord, period_px, L, meanL, N)
% BUILD_SAWTOOTH  Map projected coordinates to sawtooth luminance levels.
    pos   = mod(coord / period_px, 1);
    frame = meanL * ones(N, N);
    frame(pos < 0.25)               = L(1);
    frame(pos >= 0.25 & pos < 0.50) = L(2);
    frame(pos >= 0.50 & pos < 0.75) = L(3);
    frame(pos >= 0.75)              = L(4);
    patch = single(frame);
end
