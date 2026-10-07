function mI = make_DH_stimulus_1D(N, n_frames, framerate, dov, ...
                                   cond_idx, sigma_t_ms, phase_mov, phase_stat, ...
                                   speed_dps, sfm_list, sfs_list, meanL)
% MAKE_DH_STIMULUS_1D  Generate a Derrington & Henning stimulus as a 1D row.
%
%   mI = make_DH_stimulus_1D(N, n_frames, framerate, dov,
%                             cond_idx, sigma_t_ms, phase_mov, phase_stat,
%                             speed_dps, sfm_list, sfs_list, meanL)
%
%   Creates a compound grating with one or two sinusoidal components
%   (moving and/or static) modulated by a Gaussian temporal envelope.
%   The stimulus is 1D — identical along the vertical axis.
%
%   Conditions (selected by cond_idx):
%     1 = Simple_1mov       — 1 cpd grating, moving alone
%     2 = Simple_3mov       — 3 cpd grating, moving alone
%     3 = Complex_1mov_3stat — 1 cpd moving + 3 cpd static
%     4 = Complex_1stat_3mov — 3 cpd moving + 1 cpd static (reversal)
%     5 = Complex_1mov_3mov  — both components moving (CM)
%
%   Spatial axis: x = (0:N-1)*dov/N in degrees.
%   Temporal envelope: Gaussian centred at T/2 with sigma = sigma_t_ms.
%
% Inputs
%   N (int) number of spatial samples (pixels along x)
%   n_frames (int) number of temporal frames
%   framerate (double) frames per second [fps]
%   dov (double) field of view [degrees]
%   cond_idx (int) condition index (1–5)
%   sigma_t_ms (double) gaussian envelope width [ms]
%   phase_mov (double) initial phase of the moving component [radians]
%   phase_stat (double) initial phase of the static component [radians]
%   speed_dps (double) drift speed [degrees/second]
%   sfm_list (1x5 double) moving component SF per condition [cpd]
%   sfs_list (1x5 double) static component SF per condition [cpd]
%                             0 = no static component.
%   meanL (double) mean luminance [a.u.]
%
% Outputs
%   mI (double, n_frames x N) 1D stimulus row repeated across time
%       Luminance values in [0, 2*meanL].

    has_static = (cond_idx == 3 || cond_idx == 4);
    is_CM      = (cond_idx == 5);

    x_deg  = (0:N-1) * dov / N;
    t_vec  = 0 : n_frames-1;
    tf_fps = speed_dps / 1000;

    sfm_cpd = sfm_list(cond_idx);
    sfs_cpd = sfs_list(cond_idx);

    sigma_fr = sigma_t_ms / 1000 * framerate;
    G_env    = exp(-(t_vec - n_frames/2).^2 / (2 * sigma_fr^2));

    mI = zeros(n_frames, N);
    for tt = 1:n_frames
        row = sin(2*pi * (-sfm_cpd * x_deg + tf_fps * t_vec(tt)) + phase_mov);
        if is_CM
            row = row + sin(2*pi * (-sfs_cpd * x_deg + tf_fps * t_vec(tt)) + phase_stat);
        elseif has_static
            row = row + sin(2*pi * (-sfs_cpd * x_deg) + phase_stat);
        end
        mI(tt,:) = (1 + row * G_env(tt)) * meanL;
    end
end
