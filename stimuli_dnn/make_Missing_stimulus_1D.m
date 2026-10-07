function mI = make_Missing_stimulus_1D(N, n_frames, framerate, dov, ...
                                   sigma_t_ms, phase_mov, ...
                                   speed_dps, sfm_list, meanL, Type)
% MAKE_MISSING_STIMULUS_1D  Generate 1D Missing Fundamental 

%   mI = make_Missing_stimulus_1D(N, n_frames, framerate, dov,
%                             sigma_t_ms, phase_mov,
%                             speed_dps, sfm_list, meanL, Type)
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


  x_deg  = (0:N-1) * dov / N;
    sfm_cpd = sfm_list; 
    
    t_vec = 0 : n_frames-1;
    sigma_fr = sigma_t_ms / 1000 * framerate;
    G_env    = exp(-(t_vec - n_frames/2).^2 / (2 * sigma_fr^2));
    
    mI = zeros(n_frames, N);
    harmonics = 1:2:39;
    
    frames_per_step = round(framerate*40/1000); 
    
    dt = 1 / framerate; 
    base_dt = frames_per_step * dt; 
    
    t_stepped = 0;
    counter = 0;
    
    for tt = 1:n_frames
        counter = counter + 1;
        if counter == frames_per_step
            t_stepped = t_stepped + base_dt; 
            counter = 0;
        end
        
        row = zeros(1, N);
        
        for h = harmonics
            if Type == -1 && h == 1
                continue;
            end
           
            row = row + (4/pi)*(1/h)*sin(2*pi*h*sfm_cpd*(x_deg-speed_dps*t_stepped) + h*phase_mov);
        end
        mI(tt, :) = meanL*(1 + row*G_env(tt));
    end
end