function mI = make_Contrast_stimulus_1D(N, n_frames, framerate, dov, ...
                                   sigma_t_ms, phase_mov, ...
                                   speed_dps, sfm_list, meanL)
% MAKE_CONTRAST_STIMULUS_1D  Generate a 1D grating
%
%   mI = make_Contrast_stimulus_1D(N, n_frames, framerate, dov,
%                             sigma_t_ms, phase_mov,
%                             speed_dps, sfm_list, meanL)
%
%
%   Conditions (selected by cond_idx):
%     1 = Simple_mov       — grating moving alone

%
%   Spatial axis: x = (0:N-1)*dov/N in degrees.
%   Temporal envelope: Gaussian centred at T/2 with sigma = sigma_t_ms.
%
% Inputs
%   N (int) number of spatial samples (pixels along x)
%   n_frames (int) number of temporal frames
%   framerate (double) frames per second [fps]
%   dov (double) field of view [degrees]
%   sigma_t_ms (double) gaussian envelope width [ms]
%   phase_mov (double) initial phase of the moving component [radians]
%   speed_dps (double) drift speed [degrees/second]
%   sfm_list (1x5 double) moving component SF per condition [cpd]
%   meanL (double) mean luminance [a.u.]
%
% Outputs
%   mI (double, n_frames x N) 1D stimulus row repeated across time
%       Luminance values in [0, 2*meanL].


    x_deg  = (0:N-1) * dov / N;
    t_vec  = 0 : n_frames-1;
    sfm_cpd = sfm_list;
    sigma_fr = sigma_t_ms / 1000 * framerate;
    G_env    = exp(-(t_vec - n_frames/2).^2 / (2 * sigma_fr^2));
    mI = zeros(n_frames, N);
    for tt = 1:n_frames
        row = sin(2*pi*sfm_cpd*(x_deg - speed_dps*(t_vec(tt)/framerate)) + phase_mov);
        mI(tt,:) = (1 + row*G_env(tt))*meanL;
    end
end
