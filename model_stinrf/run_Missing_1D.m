%% run_Missing_1D.m
%
% Runs the missing fundamental experiment
%
% Each condition is run once with the spatial phase fixed at 0, the same
% phase used for the optical-flow networks, so that both see identical
% stimuli. The model response is phase-invariant here (under 0.2% change
% across phases), so a single fixed phase loses nothing and is reproducible.
% The output is the mean model response per condition (signed, rightward
% positive).
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
%   n_trials     (int)    repetitions per condition (1: phase is fixed)
%   rng_seed     (int)    random seed, 0 for unseeded
%   Type         (int)    1 square wave -1 missing fundamental
%
clear; close all;
cfg   = stINRF3_config('Cont');
model = cfg.model;

%% STIMULUS PARAMETERS
N            = 400;          % spatial samples [px]
dov          = 4.0;          % field of view [dva]
framerate    = 120;         % temporal sampling [fps] (1 ms resolution)
duration_ms  = 500;          % stimulus duration [ms]
speed_dps    = 6;            % grating drift speed [degrees/second]
meanL        = 100;          % mean luminance [a.u.]
sfm_list     = [1];          % Spatial frequency [cpd]
Type         = [1, -1];      % 1 square wave, -1 missing fundamental
sigma_t      = 10000;         % Gaussian temporal envelope widths to sweep [ms]

%% TRIAL 
n_trials = 1;               % one run, phase fixed at 0 (as for the DNNs)
rng_seed = 0;               % 0 = no seed (different each run)

%% READOUT
% Spatial readout filter c_s: central row of the 2D isotropic Gaussian,
% s_min = s_max.  sigma_c = L/sqrt(12) with L = 200 px the width of the
% sensor array it replaces.
sigma_c      = 57.74;        % [px]
make_videos  = false;        % true = also write the stimulus videos

%% DERIVED
n_frames   = round(duration_ms / 1000 * framerate);
n_type     = numel(Type);
j0         = round(N/2);
j_sensors  = (j0-99):(j0+100);  % 200 sensors centred on the stimulus
if rng_seed > 0, rng(rng_seed); end

% Pre-calculate physical units for axes labeling
x_deg = (0:N-1) * dov / N - dov/2;
t_sec = (0:n_frames-1) / framerate;

fprintf('\n==== RESULTS MISSING FUNDAMENTAL (1D) ====\n');
fprintf('Readout: 200 sensors\n');

%% RUN
fprintf('  %-12s %-16s %-12s\n', 'Type', 'Response', 'Trials');
fprintf('  %s\n', repmat('-', 1, 44));

figure(1)
set(gcf, 'Position', [50, 150, 1400, 500]);

for si = 1:n_type
    n_correct = 0;
    for tr = 1:n_trials
        phase_mov  = 0;              % fixed phase, as for the DNNs
        mI = make_Missing_stimulus_1D(N, n_frames, framerate, dov, ...
            sigma_t, phase_mov, ...
            speed_dps, sfm_list, meanL, Type(si));
        
        % Normalize and apply contrast scaling
        %mI = meanL * (1 + 2 * (mI ./ max(abs(mI(:))) - 0.5));
        
      
        subplot(1, 3, si) 
        imagesc(x_deg, t_sec, mI)
        colormap(gray)
        axis square
        xlabel('Space (degrees)')
        ylabel('Time (seconds)')
        
        if tr == n_trials && make_videos
        if Type(si) == 1
            title('Square Wave (Type = 1)')
            video_filename = 'Square_Wave_Stimulus.avi';
        else
            title('Missing Fundamental (Type = -1)')
            video_filename = 'Missing_Fundamental_Stimulus.avi';
        end

        pause(0.01)

        %%%Make a video of each type.

        v = VideoWriter(video_filename, 'Motion JPEG AVI');
        v.FrameRate = 60;   % Hz
        v.Quality = 100;   
        open(v);
        
        height_2d = 400; %%Height of the image   
        
        mI_normalized = mI - min(mI(:));
        mI_uint8 = uint8(255 * (mI_normalized / max(mI_normalized(:))));
        
        for tt = 1:n_frames
           
            frame_2d = repmat(mI_uint8(tt, :), height_2d, 1);
            writeVideo(v, frame_2d);
        end
        close(v);
        fprintf('Video saved: %s\n', video_filename);
        end

        % Compute model response
        O_sensors = stINRF3_compute_1D(mI, dov, model, framerate, N, j_sensors);
        Response(tr) = stINRF3_readout_gaussian_1D(O_sensors, j_sensors, sigma_c);
    end
    All_trials(si,:) = Response(:).';
    Mean_response(si) = mean(Response);
    SD_response(si) = std(Response);

    fprintf('  %-12s %-16s\n', ...
        sprintf('%.3f', Type(si)), ...
        sprintf('%d    %d', mean(Response), n_trials));
end
fprintf('\n');

SEM_response = SD_response ./ sqrt(n_trials); 

 subplot(1, 3, 3)
b = bar(Mean_response, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', [0 0 0], 'LineWidth', 1.5);
hold on;
errorbar(1:n_type, Mean_response, SEM_response, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 15);
xlim([0.5, n_type + 0.5]);
ax = gca;
ax.XTick = 1:n_type;
ax.XTickLabel = {'Square wave', 'Missing fundamental'};
ylabel('MS-INRF mean response (a.u.)');
title('Model Responses');
grid on;
axis square;
hold off;

%% SAVE  (small variables only, so the plotting script can reuse them)
res_dir = stinrf_results_dir(mfilename);
if isempty(res_dir), res_dir = pwd; end
if ~exist(res_dir,'dir'), mkdir(res_dir); end
wv = whos; small = {wv([wv.bytes] < 2e5).name};
save(fullfile(res_dir, ['Missing_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), small{:});
fprintf('Saved to %s\n', res_dir);

export_Missing;
