function mI = make_Interaction_stimulus_1D(N, n_frames, framerate, dov, ...
                                   cond_idx, sigma_t_ms, phase_High, phase_Low, ...
                                   v_High_dps, v_Low_dps, sfm_High_cpd, sfm_Low_cpd, ...
                                   meanL, direction)
% MAKE_INTERACTION_STIMULUS_1D  Cross-scale interaction stimulus, 1D.
%
%   Sum of a low and a high spatial frequency component under a Gaussian
%   temporal envelope. Which of the two drifts is set by cond_idx:
%     1 S1  — low frequency alone, moving
%     2 S3  — high frequency alone, moving
%     3 CS  — low moving  + high static
%     4 CSr — high moving + low static   (direction reversal expected)
%     5 CM  — both moving
%
% Inputs
%   N (int) number of spatial samples (pixels along x)
%   n_frames (int) number of temporal frames
%   framerate (double) frames per second [fps]
%   dov (double) field of view [degrees]
%   cond_idx (int) condition index 1-5
%   sigma_t_ms (double) Gaussian temporal envelope width [ms]
%   phase_High (double) initial phase of the high-frequency component [rad]
%   phase_Low (double) initial phase of the low-frequency component [rad]
%   v_High_dps (double) drift speed of the high-frequency component [deg/s]
%   v_Low_dps (double) drift speed of the low-frequency component [deg/s]
%   sfm_High_cpd (double) high component spatial frequency [cpd]
%   sfm_Low_cpd (double) low component spatial frequency [cpd]
%   meanL (double) mean luminance [a.u.]
%   direction (int, optional) +1 rightward (default), -1 leftward
%
% Outputs
%   mI (double, n_frames x N) space-time luminance matrix, rows are frames

    if nargin < 13 || isempty(direction), direction = 1; end

    move_low  = ismember(cond_idx, [1 3 5]);
    move_high = ismember(cond_idx, [2 4 5]);
    use_low   = (cond_idx ~= 2);
    use_high  = (cond_idx ~= 1);

    x_deg    = (0:N-1) * dov / N;
    t_vec    = 0 : n_frames-1;
    sigma_fr = sigma_t_ms / 1000 * framerate;
    G_env    = exp(-(t_vec - n_frames/2).^2 / (2 * sigma_fr^2));

    mI = zeros(n_frames, N);
    for tt = 1:n_frames
        tsec = t_vec(tt) / framerate;

        Stimulus = zeros(1, N);
        if use_low
            x_low = x_deg - move_low * direction * v_Low_dps * tsec;
            Stimulus = Stimulus + sin(2*pi * sfm_Low_cpd * x_low + phase_Low);
        end
        if use_high
            x_high = x_deg - move_high * direction * v_High_dps * tsec;
            Stimulus = Stimulus + sin(2*pi * sfm_High_cpd * x_high + phase_High);
        end

        mI(tt,:) = meanL * (1 + 0.5 * G_env(tt) * Stimulus);
    end
end
