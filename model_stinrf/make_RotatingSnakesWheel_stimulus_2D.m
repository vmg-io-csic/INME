function [TCB, I_lum] = make_RotatingSnakesWheel_stimulus_2D(variant, N, ...
                                    n_frames, onset_f, offset_f, meanL)
%MAKE_ROTATINGSNAKESWHEEL_STIMULUS_2D  Static snakes wheel with abrupt onset.
%
%   Builds a wheel of n_rings concentric annuli, each divided into n_units
%   identical angular units carrying the same luminance profile, with
%   successive rings rotated relative to one another. The pattern never moves:
%   it appears abruptly at frame onset_f over a uniform background and
%   disappears after frame offset_f.
%
%   For every pixel the generator computes its radius as a fraction of the
%   wheel radius and its polar angle, finds which annulus it belongs to, maps
%   the angle to a phase within one angular unit (offsetting it by the ring
%   index) and reads the luminance profile at that phase by linear
%   interpolation, wrapping around so the profile is periodic.
%
%   Inputs
%     variant   char. 'illusion' for the asymmetric wheel that elicits illusory
%               rotation, 'symmetric' for the control wheel. See
%               rotating_snakes_params for the geometry and luminance of each.
%     N         int. Canvas side in pixels. The canvas is square and the wheel
%               fills it, with a radius of 0.49*N.
%     n_frames  int. Total number of frames in the sequence.
%     onset_f   int. First frame in which the wheel is visible.
%     offset_f  int. Last frame in which the wheel is visible.
%     meanL     double. Background luminance, in the same arbitrary units as
%               the profile, which spans [0, 100].
%
%   Outputs
%     TCB    single, N x N x n_frames. Stimulus tensor, ready for
%            stINRF3_compute.
%     I_lum  double, N x N. The static pattern on its own, in [0, 100] inside
%            the wheel and meanL outside it.
%
%   Example
%     [TCB, I] = make_RotatingSnakesWheel_stimulus_2D('illusion', 512, 60, 11, 60, 50);

    p = rotating_snakes_params(variant);

    cx = (N+1)/2;
    cy = (N+1)/2;
    r_out = 0.49 * N;             % the wheel fills the canvas

    [Xg, Yg] = meshgrid(1:N, 1:N);
    dx = Xg - cx;
    dy = Yg - cy;
    r_frac = sqrt(dx.^2 + dy.^2) / r_out;
    theta  = mod(atan2(dy, dx), 2*pi);

    ring_idx = discretize(r_frac, p.edges);   % NaN outside the wheel

    n_prof   = numel(p.unit_profile);
    prof_ext = [p.unit_profile, p.unit_profile(1)];   % wrap for periodic interp
    u_grid   = linspace(0, 1, n_prof + 1);

    I_lum = meanL * ones(N, N);
    valid = ~isnan(ring_idx) & (r_frac <= 1);
    k = ring_idx(valid);
    u = mod(theta(valid)/(2*pi) * p.n_units - (k - 1) * p.angular_offset, 1);
    I_lum(valid) = interp1(u_grid, prof_ext, u, 'linear');

    TCB   = single(meanL * ones(N, N, n_frames));
    patch = single(I_lum);
    for t = onset_f : min(offset_f, n_frames)
        TCB(:,:,t) = patch;
    end
end
