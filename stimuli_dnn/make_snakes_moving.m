function TCB = make_snakes_moving(N, n_frames, onset_f, offset_f, ...
                                   period_px, orientation_deg, L, meanL, ...
                                   vel_dps, dpi, framerate, phase_shift_px)
% MAKE_SNAKES_MOVING  Generate a physically drifting sawtooth pattern.
%
%   TCB = make_snakes_moving(N, n_frames, onset_f, offset_f,
%                             period_px, orientation_deg, L, meanL,
%                             vel_dps, dpi, framerate, phase_shift_px)
%
%   Same spatial structure as make_snakes_static, but the pattern
%   translates at vel_dps degrees/second from frame to frame.  Used as
%   a real-motion control condition to verify that the model responds
%   to genuine motion (as opposed to illusory motion from luminance
%   asymmetry at onset).
%
% Inputs
%   N (int) image size [px] (square: N x N)
%   n_frames (int) total number of frames
%   onset_f (int) first visible frame
%   offset_f (int) last visible frame
%   period_px (double) spatial period of the sawtooth [px]
%   orientation_deg (double) pattern orientation [degrees]
%   L (1x4 double) luminance levels per quarter-period
%   meanL (double) background luminance [a.u.]
%   vel_dps (double) drift velocity [degrees/second]
%   dpi (double) pixels per degree (= N / dov)
%   framerate (double) frames per second [fps]
%   phase_shift_px (double) initial spatial phase offset [px] (default: 0)
%
% Outputs
%   TCB (single, N x N x n_frames) stimulus tensor

    if nargin < 12, phase_shift_px = 0; end

    vel_ppf = vel_dps * dpi / framerate;

    [X, Y] = meshgrid(0:N-1, 0:N-1);
    X = X - (N-1)/2;
    Y = Y - (N-1)/2;
    coord_base = -X*cosd(orientation_deg) + Y*sind(orientation_deg);

    TCB = single(meanL * ones(N, N, n_frames));
    for t = onset_f : min(offset_f, n_frames)
        shift = phase_shift_px + (t - onset_f) * vel_ppf;
        TCB(:,:,t) = build_sawtooth(coord_base + shift, period_px, L, meanL, N);
    end
end

function patch = build_sawtooth(coord, period_px, L, meanL, N)
% BUILD_SAWTOOTH  Map projected coordinates to sawtooth luminance levels.
%   See make_snakes_static.m > build_sawtooth for documentation.
    pos   = mod(coord / period_px, 1);
    frame = meanL * ones(N, N);
    frame(pos < 0.25)               = L(1);
    frame(pos >= 0.25 & pos < 0.50) = L(2);
    frame(pos >= 0.50 & pos < 0.75) = L(3);
    frame(pos >= 0.75)              = L(4);
    patch = single(frame);
end
