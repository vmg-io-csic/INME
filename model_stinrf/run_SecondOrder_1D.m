%% run_SecondOrder_1D.m
%
% Runs the second order experiment
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
%   Type         (int)    1 square wave -1 missing fundamental
%
clear; close all;
cfg   = stINRF3_config('Cont');
model = cfg.model;
% Second-order motion only, as run_LightEdges.m and run_Uneven.m set their
% own temporal parameters. With the shared sigmaMu = 0.05 deg the A-stage
% smoothing wipes out most of the 1 px carrier before the nonlinearity, so
% there is little texture left to rectify and the response to each trial is
% dominated by the first-order energy of that particular carrier (64% of
% trials in the right direction). With these three values every trial comes
% out in the direction of the envelope.
model.sigmaMu = 0.02;        % [dva] (shared set: 0.05)
model.p       = 0.1;         % []    (shared set: 0.2)
model.Nm      = 60;          % [ms]  (shared set: 125)

%% STIMULUS PARAMETERS
N            = 400;          % spatial samples [px]
dov          = 4.0;          % field of view [dva]
framerate    = 120;         % temporal sampling [fps] (1 ms resolution)
duration_ms  = 500;          % stimulus duration [ms]
speed_dps    = 8;            % grating drift speed [degrees/second]
meanL        = 100;          % mean luminance [a.u.]
sfm_envelope_cpd  = 1;       % Envelope Spatial frequency [cpd]
sfm_carrier_cpd   = 4;       % Carrier Spatial frequency [cpd] (for Carrier_type=0)
Carrier_type = 1;            % 1 1D static noise  0 grating. The figure uses the noise carrier: a static sinusoidal carrier adds two first-order components drifting in opposite directions (fc +/- fe)
Type         = [-1, 1];      % -1 moving left, 1 moving right
sigma_t      = 10000;         % Gaussian temporal envelope widths to sweep [ms]

%% TRIAL 
n_trials = 100;               % random-phase repetitions
rng_seed = 1;               % 0 = no seed (different each run); 1 gives the figure's carriers and phases

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

fprintf('\n==== RESULTS SECOND-ORDER MOTION (1D) ====\n');
fprintf('Readout: 200 sensors\n');

%% RUN
fprintf('  %-12s %-16s %-12s\n', 'Type', 'Response', 'Trials');
fprintf('  %s\n', repmat('-', 1, 44));


figure(1)
set(gcf, 'Position', [50, 150, 1400, 500]);

for si = 1:n_type
    n_correct = 0;
    for tr = 1:n_trials
        phase_mov  = rand * 2*pi;
        phase_carrier  = rand * 2*pi;

        mI = make_SecondOrder_stimulus_1D(N, n_frames, framerate, dov, ...
            sigma_t, phase_mov, phase_carrier, ...
            speed_dps, sfm_envelope_cpd, sfm_carrier_cpd, meanL, Type(si), Carrier_type);
        
           
        subplot(1, 3, si) 
        imagesc(x_deg, t_sec, mI)
            if Type(si) == -1
            title('Second Order (Left) (Type = -1)')
            else
            title('Second Order (Right) (Type = 1)')
            end
        colormap(gray)
        axis square
        xlabel('Space (degrees)')
        ylabel('Time (seconds)')
        pause(0.01)
        
        if tr == n_trials && make_videos
        if Type(si) == -1
            video_filename = 'Second_order_left.avi';
        else
            video_filename = 'Second_order_right.avi';
        end



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
        sprintf('%d      %d', mean(Response), n_trials));
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
ax.XTickLabel = {'Left', 'Right'};
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
save(fullfile(res_dir, ['SecondOrder_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), small{:});
fprintf('Saved to %s\n', res_dir);

export_SecondOrder;
