function mI = make_Opponency_stimulus_1D(N, n_frames, framerate, dov, ...
            sigma_t_ms, speed_dps, bar_width, meanL, direction, inver)
% MAKE_OPPONENCY_STIMULUS_1D  Generate a 1D drifting white bar on a black background
%
% Inputs
%   N (int) number of spatial samples (pixels along x)
%   n_frames (int) number of temporal frames
%   framerate (double) frames per second [fps]
%   dov (double) field of view [degrees]
%   sigma_t (double) gaussian envelope width [seconds]
%   speed_dps (double) drift speed [degrees/second]
%   bar_width (double) width of the bar [degrees]
%   meanL (double) maximum luminance of the white bar [a.u.]
%   direction (int)  1 for moving right, -1 for moving left
%   inver  (int) 0 (white bar, black background) 1 (black bar, white background)
% Outputs
%   mI (double, n_frames x N) 1D stimulus matrix across time


    x_deg = (0:N-1) * dov / N;
    t_vec = (0:n_frames-1); 
    
    t_center = n_frames / 2;
    sigma_fr = sigma_t_ms / 1000 * framerate;
    G_env = exp(-(t_vec - t_center).^2 / (2 * sigma_fr^2));
    mI = zeros(n_frames, N);
    
    % Always calculate moving from left (0) to right
    start_pos = 0; 
    
    for tt = 1:n_frames
        % Linear motion forward
        bar_center = start_pos + speed_dps * t_vec(tt)/1000;
        
        % Create the bar mask
        bar_mask = abs(x_deg - bar_center) <= (bar_width / 2);
        
        mI(tt, bar_mask) = meanL * G_env(tt);
    end
    
    %If direction is -1, flip the space dimension (columns) horizontally
    if direction == -1
        mI = fliplr(mI);
    end
if inver == 1
    mI = meanL-mI;
end

end