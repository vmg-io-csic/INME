function mI = make_masking_stimulus_1D(N, n_frames, framerate, dov, ...
                                   sigma_t_ms, phase_mov, ...
                                   speed_dps, sfm_list,sfmask_cpd, meanL, masking, Contrast_signal, Contrast_mask, temporal_frequency_mask)
% MAKE_MASKING_STIMULUS_1D  Generate a 1D grating
%
%   mI = make_Masking_stimulus_1D(N, n_frames, framerate, dov,
%                             sigma_t_ms, phase_mov, phase_stat,
%                             speed_dps, sfm_list, meanL, mask,Contrast_signal, Contrast_mask)
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
        signal = Contrast_signal * sin(2*pi*sfm_cpd*(x_deg - speed_dps * (t_vec(tt)/framerate)) + phase_mov);
        mask   = Contrast_mask * sin(2*pi*(sfmask_cpd * x_deg + (sign(rand - 0.5))*temporal_frequency_mask * (t_vec(tt)/framerate)));
        
        if masking == 1
            row = (signal + mask)/2;
        else
            row = signal;
        end
        mI(tt,:) = (1 + row * G_env(tt)) * meanL;
    end
end
