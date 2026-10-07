%% run_Masking_1D.m
%
% Runs the masking experiment
%
% For each sigma_t, n_trials are run with random spatial phases.
% The output is the percentage of trials where the model reports
% the correct (positive) direction.
%
% Inputs
%   N            (int)    number of spatial pixels
%   dov          (double) field of view [dva]
%   framerate    (int)    frame rate [fps]
%   duration_ms  (int)    total duration [ms]
%   speed_dps    (double) grating drift speed in degrees/second
%   meanL        (double) mean luminance [a.u.]
%   sfm_list     (1)      moving component SF per condition in cpd
%   sigma_t_list (double) Gaussian envelope width to sweep in ms
%   sfmask_cpd   (double) Spatial frequency of the mask, cpd
%   condition    (int)    condition index 1-5 (see CONDITION SELECTOR)
%   n_trials     (int)    random-phase repetitions per sigma_t
%   rng_seed     (int)    random seed, 0 for unseeded
%   Contrast     (int)    Contrast of the grating
%

clear; close all;

cfg   = stINRF3_config('Cont');
% The model keeps its default sf0 and the signal is at 1 cpd, as in
% Programs_ISP. An earlier version tuned the detector to a 2.5 cpd signal
% (cfg.model.sf0 = 2.5).
model = cfg.model;

%% STIMULUS PARAMETERS
N            = 400;          % spatial samples [px]
dov          = 4.0;          % field of view [dva]
framerate    = 120;         % temporal sampling [fps] (1 ms resolution)
duration_ms  = 500;          % stimulus duration [ms]
speed_dps    = 10;           % grating drift speed [degrees/second]
meanL        = 100;          % mean luminance [a.u.]
temporal_frequency_mask  = 10; %Hz

sfm_list = 1; %Signal Spatial frequency [cpd] (as in Programs_ISP)
sfmask_cpd = [0.02 0.05 0.2 0.25 0.5 1 1.5 2 2.5 3 4 6 12];%%% Mask Spatial frequency [cpd] (Programs_ISP set plus 0.25 1.5 2 2.5 4)

masking = 1; % 1 Yes  0 No masking

Response_no_mask= 5.636215e+00; %this is the response to simple grating without mask (framerate=120). If you change parameters you have to compute this value again

Contrast_signal = 0.8;
Contrast_mask = 0.8;

% Gaussian temporal envelope widths to sweep [ms]
sigma_t=1000;


%% TRIAL / READOUT
n_trials = 50;               % random-phase repetitions
rng_seed = 0;                % 0 = no seed (different each run)

%% READOUT
% Spatial readout filter c_s: central row of the 2D isotropic Gaussian,
% s_min = s_max.  sigma_c = L/sqrt(12) with L = 200 px the width of the
% sensor array it replaces.
sigma_c      = 57.74;        % [px]
make_videos  = false;        % true = also write the stimulus videos

%% DERIVED
n_frames   = round(duration_ms / 1000 * framerate);
n_masking   = numel(sfmask_cpd);
j0         = round(N/2);
j_sensors  = (j0-99):(j0+100);  % 200 sensors centred on the stimulus

if rng_seed > 0, rng(rng_seed); end

fprintf('\n==== RESULTS MASKING EXPERIMENT (1D) ====\n');
fprintf('Readout: 200 sensors\n');

%% RUN
fprintf('  %-12s %-16s %-12s\n', 'Spatial freq. mask', 'Response', 'Trials');
fprintf('  %s\n', repmat('-', 1, 44));

% --- Calculate Physical Coordinate Vectors ---
% Space and Time coordinates
x_deg = (0:N-1) * dov / N - dov/2;                      % Degrees of visual angle [dva]
t_sec = (0:n_frames-1) / framerate;             % Time [seconds]

% Frequency coordinates (Fourier domain)
dx = dov / N;                                   % Spatial sampling interval [degrees]
sf_nyquist = 1 / (2 * dx);                      % Nyquist spatial frequency [cpd]
fx_cpd = linspace(-sf_nyquist, sf_nyquist, N);  % Spatial frequency vector

dt = 1 / framerate;                             % Temporal sampling interval [seconds]
tf_nyquist = 1 / (2 * dt);                      % Nyquist temporal frequency [Hz]
ft_hz = linspace(-tf_nyquist, tf_nyquist, n_frames); % Temporal frequency vector

for si = 1:n_masking
    n_correct = 0;
    for tr = 1:n_trials
        phase_mov  = rand * 2*pi;
        mI = make_Masking_stimulus_1D(N, n_frames, framerate, dov, ...
            sigma_t, phase_mov, ...
            speed_dps, sfm_list, sfmask_cpd(si), meanL, masking, Contrast_signal, Contrast_mask, temporal_frequency_mask);
        
       
       % mI = meanL * (1 + Contrast * 2 * (mI ./ max(abs(mI(:))) - 0.5));
        
        %2D Fourier Amplitude Spectrum
        mI_fft = fftshift(fft2(mI - mean(mI(:)))); 
        mI_amp = abs(mI_fft); 
        

        figure(1)
        set(gcf, 'Position', [100, 200, 1000, 450]); 
        
        %Space-Time Image (mI)
        subplot(1, 2, 1)
        imagesc(x_deg, t_sec, mI)
        colormap(gca, gray)
        title(sprintf('Stimulus mI (sf_{mask} = %.2f cpd)', sfmask_cpd(si)))
        xlabel('Space (degrees)')
        ylabel('Time (seconds)')
        axis square
        
        % Fourier Amplitude Spectrum
        subplot(1, 2, 2)
        imagesc(fx_cpd, ft_hz, mI_amp) 
        colormap(gca, jet) 
        title(sprintf('Fourier Spectrum (sf_{mask} = %.2f cpd)', sfmask_cpd(si)))
        xlabel('Spatial Frequency (cpd)')
        ylabel('Temporal Frequency (Hz)')
        axis square
        set(gca, 'YDir', 'normal')
        ylim([-40 40]); 
        xlim([-15 15]);   
        
        pause(0.01)
        
        % Model response
        O_sensors = stINRF3_compute_1D(mI, dov, model, framerate, N, j_sensors);
        Response(tr) = stINRF3_readout_gaussian_1D(O_sensors, j_sensors, sigma_c);
    end
    Mean_response(si) = mean(Response);
    fprintf('  %-12s %-16s\n', ...
        sprintf('%.3f', sfmask_cpd(si)), ...
        sprintf('%d          %d', mean(Response), n_trials));
end

fprintf('\n');

figure(3)
 set(gcf, 'Position', [200, 250, 1000, 450]); 
subplot(1,2,1)
semilogx(sfmask_cpd, Mean_response, 'bo-', [0.01, 32], [Response_no_mask, Response_no_mask], 'r--', 'LineWidth', 2)
title('Plot linear-log')
ax = gca; 
ax.XTick = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 4 8 16 32]; 
ax.XTickLabel = {'0.01','0.02','0.05','0.1','0.2', '0.5', '1', '2', '4', '8', '16', '32'}; 
ax.YTick = [0, 1, 2, 3, 4, 5, 6, 7]; 
ax.YTickLabel = {'0', '1', '2', '3', '4', '5', '6', '7'};
ylim([0 7]); 
xlim([0.01 32]); 
ylabel('MS-INRF response (a.u.)')
xlabel('Spatial frequency of the mask (cpd)')
axis square
grid on

subplot(1,2,2)
plot(sfmask_cpd, Mean_response, 'bo-', [0.01, 32], [Response_no_mask, Response_no_mask], 'r--', 'LineWidth', 2)
title('Plot linear-linear') 
ax = gca; 
ax.XTick = [0.5 1 2 4 8 12];
ax.XTickLabel = {'0.5', '1', '2', '4', '8', '12'}; 
ax.YTick = [0, 1, 2, 3, 4, 5, 6, 7]; 
ax.YTickLabel = {'0', '1', '2', '3', '4', '5', '6', '7'};
ylim([0 7]); 
xlim([0.01 12]); 
ylabel('MS-INRF response (a.u.)')
xlabel('Spatial frequency of the mask (cpd)')
axis square
grid on

%% SAVE  (small variables only, so the plotting script can reuse them)
res_dir = stinrf_results_dir(mfilename);
if isempty(res_dir), res_dir = pwd; end
if ~exist(res_dir,'dir'), mkdir(res_dir); end
wv = whos; small = {wv([wv.bytes] < 2e5).name};
save(fullfile(res_dir, ['Masking_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), small{:});
fprintf('Saved to %s\n', res_dir);

export_Masking;
