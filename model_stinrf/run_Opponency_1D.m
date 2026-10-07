%% run_Opponency_1D.m
%
% Runs the opponency experiment
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
%   condition    (int)    condition index 1-5 (see CONDITION SELECTOR)
%   n_trials     (int)    random-phase repetitions per sigma_t
%   rng_seed     (int)    random seed, 0 for unseeded
%   Contrast     (int)    Contrast of the grating
%   inver  (int) 1 (white bar, black background) 0 (black bar, white background)
%

clear; close all;

cfg   = stINRF3_config('Cont');
model = cfg.model;

%% STIMULUS PARAMETERS
N            = 400;          % spatial samples [px]
dov          = 4.0;          % field of view [dva]
framerate    = 1000;         % temporal sampling [fps] (1 ms resolution)
duration_ms  = 500;          % stimulus duration [ms]
speed_dps    = 8;           % grating drift speed [degrees/second]
meanL        = 100;          % mean luminance [a.u.]

% Bar width [degrees]
bar_width = [0.25]; 


% Gaussian temporal envelope widths to sweep [ms]
sigma_t=100000;


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
j0         = round(N/2);
j_sensors  = (j0-99):(j0+100);  % 200 sensors centred on the stimulus

if rng_seed > 0, rng(rng_seed); end

% Pre-calculate physical units for axes labeling
x_deg = (0:N-1) * dov / N - dov/2;
t_sec = (0:n_frames-1) / framerate;


fprintf('\n==== RESULTS OPPONENCY (1D) ====\n');
fprintf('Readout: 200 sensors\n');

%% RUN
fprintf('  %-12s %-16s %-12s\n', 'Direction', 'Response', 'Trials');
fprintf('  %s\n', repmat('-', 1, 44));

Direction=[1 -1 1 -1];

figure(1)
set(gcf, 'Position', [50, 250, 1400, 350]); 

for si = 1:4
    if si > 2
        inver = 1;  %% 1 (white bar, black background)
    else 
        inver = 0;  %% 0 (black bar, white background)
    end 
    
    for tr = 1:n_trials
        mI = make_Opponency_stimulus_1D(N, n_frames, framerate, dov, ...
            sigma_t, speed_dps, bar_width, meanL, Direction(si), inver);
        mI = meanL*(1 + 2*(mI./abs(max(max(mI))) - 0.5));
        
       
        subplot(1, 5, si)       
        imagesc(x_deg, t_sec, mI)
        colormap(gray)
        axis square
        xlabel('Space (deg)')
        ylabel('Time (s)')
        
        if inver == 0
            bg_txt = 'White Bar';
        else
            bg_txt = 'Black Bar';
        end
        if Direction(si) == 1
            dir_txt = 'Right';
        else
            dir_txt = 'Left';
        end
        title(sprintf('%s (%s)', dir_txt, bg_txt));
        
        pause(0.01)
        
        O_sensors = stINRF3_compute_1D(mI, dov, model, framerate, N, j_sensors);
        All_Responses(tr) = stINRF3_readout_gaussian_1D(O_sensors, j_sensors, sigma_c);
    end
    
    All_trials(si,:) = All_Responses(:).';
    
    Mean_response(si) = mean(All_Responses);
    fprintf('  %-12s %-16s\n', ...
        sprintf('%.3f', Direction(si)), ...
        sprintf('%.3f              %d', Mean_response(si), n_trials));
end
fprintf('\n');

subplot(1, 5, 5)
b = bar(Mean_response, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', [0 0 0], 'LineWidth', 1.2);
hold on;

if n_trials > 1
    SD_response = std(All_Responses, 0, 2);
    SEM_response = SD_response ./ sqrt(n_trials);
    errorbar(1:4, Mean_response, SEM_response, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
end

ax = gca;
ax.XTick = 1:4;
ax.XTickLabel = {'R (White)', 'L (White)', 'R (Black)', 'L (Black)'};
xtickangle(45);
xlim([0.5, 4.5]); 

ylabel('MS-INRF mean response (a.u.)');
title('Model Output');
grid on;
axis square;
hold off;
%% SAVE  (small variables only, so the plotting script can reuse them)
res_dir = stinrf_results_dir(mfilename);
if isempty(res_dir), res_dir = pwd; end
if ~exist(res_dir,'dir'), mkdir(res_dir); end
wv = whos; small = {wv([wv.bytes] < 2e5).name};
save(fullfile(res_dir, ['Opponency_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), small{:});
fprintf('Saved to %s\n', res_dir);

export_Opponency;
