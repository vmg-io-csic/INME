function mI = make_SecondOrder_stimulus_1D(N, n_frames, framerate, dov, ...
                                   sigma_t_ms, phase_mov, phase_carrier, ...
                                   speed_dps, sfm_envelope_cpd, sfm_carrier_cpd, meanL, Type, Carrier_type)
% MAKE_SECONDORDER_STIMULUS_1D  Generate 1D Missing Fundamental 

%   mI = make_SecondOrder_stimulus_1D(N, n_frames, framerate, dov,
%                             sigma_t_ms, phase_mov,
%                             speed_dps, sfm_list, meanL, Type, Carrier_type)
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
%   sfm_envelope_cpd (double) moving component SF per condition [cpd]
%   sfm_carrier_cpd (double) carrier SF [cpd]
%   meanL (double) mean luminance [a.u.]
%   Type (int) -1 leftward motion, 1 for rightward
%
% Outputs
%   mI (double, n_frames x N) 1D stimulus row repeated across time
%       Luminance values in [0, 2*meanL].

    
    x_deg  = (0:N-1) * dov / N;
    t_vec  = 0 : n_frames-1;
    
    if Carrier_type == 1
        Carrier = 2 * (rand(1, N) > 0.5) - 1; %% Static carrier
    end
    
    mI = zeros(n_frames, N);
    for tt = 1:n_frames
        if Type==-1; direction=-1; else; direction=1; end
        
         Envelope = sin(2*pi*sfm_envelope_cpd*(x_deg - direction*speed_dps*(t_vec(tt)/framerate)) + phase_mov);
        
        if Carrier_type==0; Carrier = sin(2*pi*sfm_carrier_cpd*x_deg + phase_carrier); end
        
        %Carrier = 2 * (rand(1, N) > 0.5) - 1; %Dynamic carrier
        mI(tt,:) = meanL*(1 + 0.5*(1+1*Envelope).*Carrier);
        %mI(tt,:) = meanL*(1 + 0.5*Envelope); %%Test for Luminance
    end