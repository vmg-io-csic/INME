function mI = make_ReversePhi_stimulus_1D(N, n_frames, framerate, dov, ...
                                   sigma_t_ms, phase_mov, ...
                                   speed_dps, meanL, Type, Contrast)
% MAKE_REVERSEPHI_STIMULUS_1D  Generate 1D Missing Fundamental with stepped time
%   mI = make_Missing_stimulus_1D(N, n_frames, framerate, dov,
%                             sigma_t_ms, phase_mov,
%                             speed_dps, meanL, Type, Contrast)
%
%
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
%   Type (int) 1 for Square Wave, -1 for Missing Fundamental
%
% Outputs
%   mI (double, n_frames x N) 1D stimulus row repeated across time
%       Luminance values in [0, 2*meanL].


  x_deg = (0:N-1) * dov / N;
    dx = dov / N; % Tamaño del píxel en grados
    
    t_vec = 0 : n_frames-1;
    sigma_fr = sigma_t_ms / 1000 * framerate;
    G_env = exp(-(t_vec - n_frames/2).^2 / (2 * sigma_fr^2));
    
    mI = zeros(n_frames, N);

    rng('shuffle'); 
    base_noise = sign(randn(1, N)); 
    
    dt = 1 / framerate; %
    
   
    current_contrast_sign = 1;
    last_pixel_shift = 0; 
    
    for tt = 1:n_frames

        t_real = tt * dt;
        
        spatial_shift_deg = speed_dps * t_real;
        pixel_shift = round(spatial_shift_deg / dx); 

        if Type == -1
            if pixel_shift ~= last_pixel_shift
                current_contrast_sign = -current_contrast_sign; % 
                last_pixel_shift = pixel_shift; % 
            end
        else
            current_contrast_sign = 1;
        end
        
        shifted_noise = circshift(base_noise, pixel_shift);

        contrast_modulation = Contrast * shifted_noise * current_contrast_sign;
        mI(tt, :) = meanL * (1 + contrast_modulation * G_env(tt));
    end
end