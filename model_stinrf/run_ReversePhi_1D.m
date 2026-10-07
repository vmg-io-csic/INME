%% run_ReversePhi_1D.m
%
% Runs the reverse phi motion experiment
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
%   Type         (int)    1 Phi motion -1 Reverse Phi motion
%
clear; close all;
cfg   = stINRF3_config('Cont');
model = cfg.model;

%% STIMULUS PARAMETERS
N = 400;              % spatial samples [px]
dov = 16;             % field of view [dva]
framerate = 120;      % temporal sampling [fps] (1 ms resolution)
duration_ms = 500;    % stimulus duration [ms]
meanL = 100;          % mean luminance [a.u.]
Type = [1, -1];       % 1 square wave, -1 missing fundamental
sigma_t = 10000;      % Gaussian temporal envelope widths to sweep [ms]
Contrast = 0.8;       % Reducing contrast make the responses more similar for each type
P = 4;                       % Pixels per frame
speed_dps = P * (dov / N) * framerate;  % Speed in deg/s needed to move each frame 1 pixel


%% TRIAL / READOUT
n_trials = 20;               % random-phase repetitions
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

fprintf('\n==== RESULTS REVERSE PHI (1D) ====\n');
fprintf('Readout: 200 sensors\n');

%% RUN
fprintf('  %-12s %-16s %-12s\n', 'Type', 'Response', 'Trials');
fprintf('  %s\n', repmat('-', 1, 44));

figure(1)
set(gcf, 'Position', [100, 200, 950, 450]); 

for si = 1:n_type
    n_correct = 0;
    for tr = 1:n_trials
        phase_mov  = rand * 2*pi;
        mI = make_ReversePhi_stimulus_1D(N, n_frames, framerate, dov, ...
            sigma_t, phase_mov, ...
            speed_dps, meanL, Type(si), Contrast);
        
        % Normalize and apply contrast scaling
      %%%  mI2 = meanL * (1 + Contrast* 2 * (mI ./ max(abs(mI(:))) - 0.5));
        
      
        subplot(1, 2, si) 
        imagesc(x_deg, t_sec, mI)
        colormap(gray)
        axis square
        xlabel('Space (degrees)')
        ylabel('Time (seconds)')
        
        if tr == n_trials && make_videos
        if Type(si) == 1
            title('Phi motion (Type = 1)')
            video_filename = 'Phi_motion_Stimulus.avi';
        else
            title('Reverse Phi motion (Type = -1)')
            video_filename = 'Reverse_Phi_motion_Stimulus.avi';
        end



        %%%Make a video of each type.

        v = VideoWriter(video_filename, 'Motion JPEG AVI');
        v.FrameRate = 60;   % Hz
        v.Quality = 100;   
        open(v);
        
        height_2d = 400; %%Height of the image   
        
        mI_normalized = mI - min(mI(:));
        mI_uint8 = uint8(255*(mI_normalized / max(mI_normalized(:))));
        
        for tt = 1:n_frames
           
            frame_2d = repmat(mI_uint8(tt, :), height_2d, 1);
            writeVideo(v, frame_2d);
        end
        close(v);
        fprintf('Video saved: %s\n', video_filename);
        end
        pause(0.01)
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

figure(2)
b = bar(Mean_response, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', [0 0 0], 'LineWidth', 1.5);
hold on;
errorbar(1:n_type, Mean_response, SEM_response, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 15);
xlim([0.5, n_type + 0.5]);
ax = gca;
ax.XTick = 1:n_type;
ax.XTickLabel = {'Phi Motion', 'Reverse Phi'};
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
save(fullfile(res_dir, ['ReversePhi_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), small{:});
fprintf('Saved to %s\n', res_dir);

export_ReversePhi;
