function mI = make_DH_stimulus_1D_v25(N, n_frames, framerate, dov, ...
        cond_idx, sigma_t_ms, phase_mov, phase_stat, tf_hz, ...
        sfm_list, sfs_list, meanL)
% MAKE_DH_STIMULUS_1D_V25  Derrington & Henning stimulus, framerate-correct.
%
%   Identical output to make_DH_stimulus_1D.m when framerate = 1000 and
%   tf_hz = 12. The original hardcoded 1000 in
%
%       tf_fps = speed_dps / 1000;                       % <- bug B2
%
%   so the drift did not depend on `framerate` at all: the stimulus always
%   advanced 0.012 cycles/frame no matter how it was sampled, which made it
%   impossible to move the interframe displacement into the range the
%   optical-flow networks work in.
%
%   The parameter is also renamed. The original called it `speed_dps` and
%   documented it as deg/s, but the formula omits sfm_cpd, so what it sets
%   is a TEMPORAL FREQUENCY in Hz, common to both components. With
%   tf_hz = 12 the 1 cpd component drifts at 12 deg/s and the 3 cpd one at
%   4 deg/s. That is Derrington & Henning's own convention (they equate
%   temporal frequency, not speed) and it is what stINRF3 was run with, so
%   the name is corrected and the behaviour preserved.
%
% Inputs
%   N (int) spatial samples [px]
%   n_frames (int) frames
%   framerate (double) [fps]
%   dov (double) field of view [deg]
%   cond_idx (int) 1..5 -> S1 S3 CS CSr CM
%   sigma_t_ms (double) Gaussian temporal envelope sigma [ms]
%   phase_mov, phase_stat (double) initial phases [rad]
%   tf_hz (double) temporal frequency of the moving component(s) [Hz]
%   sfm_list, sfs_list (1x5 double) moving / static SF per condition [cpd]
%   meanL (double) mean luminance
%
% Outputs
%   mI (double, n_frames x N)

    has_static = (cond_idx == 3 || cond_idx == 4);
    is_CM      = (cond_idx == 5);

    x_deg  = (0:N-1) * dov / N;
    t_vec  = 0 : n_frames-1;
    tf_cpf = tf_hz / framerate;        % cycles per frame  (was tf/1000)

    sfm_cpd = sfm_list(cond_idx);
    sfs_cpd = sfs_list(cond_idx);

    sigma_fr = sigma_t_ms / 1000 * framerate;
    G_env    = exp(-(t_vec - n_frames/2).^2 / (2 * sigma_fr^2));

    mI = zeros(n_frames, N);
    for tt = 1:n_frames
        row = sin(2*pi * (-sfm_cpd * x_deg + tf_cpf * t_vec(tt)) + phase_mov);
        if is_CM
            row = row + sin(2*pi * (-sfs_cpd * x_deg + tf_cpf * t_vec(tt)) + phase_stat);
        elseif has_static
            row = row + sin(2*pi * (-sfs_cpd * x_deg) + phase_stat);
        end
        mI(tt,:) = (1 + row * G_env(tt)) * meanL;
    end
end
