%% run_TimeReversal_1D.m
%
% Runs the time-reversal symmetry-breaking experiment using the 1D
% model (stINRF3_compute_1D).
%
% Two periodic profiles drift at their paper-defined velocity/period:
%   LightEdges       — single sharp transition per period.
%   UnevenSquareWave — bars in ratio 5:5:25:25 per period.
% For each, the forward stimulus and its time-reversed version (frame
% order reversed, NOT drift direction flipped) are run through the
% model. If the model were perfectly antisymmetric under time reversal,
% R_fwd = -R_rev exactly. Any departure is quantified as:
%   V = (R_fwd - R_rev) / 2   veridical (antisymmetric) motion component
%   I = (R_fwd + R_rev) / 2   symmetry-breaking component
%   ratio = |V| / |I|
% Reference fly ratios: LightEdges = 5.72, UnevenSquareWave = 12.74.
%
% For each stimulus type, the stimulus is displayed and saved as a
% video, and the |V|/|I| ratio vs duration is plotted and saved to
% results_TimeReversal_1D/.
%
% Inputs
%   N            (int)    number of spatial pixels
%   dov          (double) field of view [dva]
%   framerate    (int)    frame rate [fps]
%   meanL        (double) mean luminance [a.u.]
%   contrast     (double) Michelson-like contrast on the ±1 profile
%   phase0_rad   (double) initial spatial phase [radians]
%   v_LightEdges, period_LightEdges (double) velocity [deg/s] and period [deg]
%   v_UnevenSquareWave, period_UnevenSquareWave (double) velocity [deg/s] and period [deg]
%   dur_ms_list  (1xD double) durations to sweep [ms]
%   ratio_fly_LightEdges, ratio_fly_UnevenSquareWave (double) reference fly |V|/|I|
%

clear; clc; close all;

cfg   = stINRF3_config('timereversal');
model = cfg.model;

%% STIMULUS PARAMETERS
N            = cfg.stim.N;
dov          = cfg.stim.dov;
framerate    = cfg.stim.framerate;
meanL        = cfg.stim.meanL;
contrast     = cfg.stim.contrast;
phase0_rad   = cfg.stim.phase0_rad;

v_LightEdges       = cfg.stim.v_LightEdges;
period_LightEdges  = cfg.stim.period_LightEdges;
v_UnevenSquareWave      = cfg.stim.v_UnevenSquareWave;
period_UnevenSquareWave = cfg.stim.period_UnevenSquareWave;
dur_ms_list  = cfg.stim.dur_ms_list;
ratio_fly_LightEdges       = cfg.stim.ratio_fly_LightEdges;
ratio_fly_UnevenSquareWave = cfg.stim.ratio_fly_UnevenSquareWave;

stim_types  = {'LightEdges', 'UnevenSquareWave'};
v_list      = [v_LightEdges, v_UnevenSquareWave];
period_list = [period_LightEdges, period_UnevenSquareWave];
ratio_fly_l = [ratio_fly_LightEdges, ratio_fly_UnevenSquareWave];


%% READOUT
% Spatial readout filter c_s: central row of the 2D isotropic Gaussian,
% s_min = s_max.  sigma_c = L/sqrt(12) with L = 200 px the width of the
% sensor array it replaces.
sigma_c      = 57.74;        % [px]
make_videos  = false;        % true = also write the stimulus videos

%% TRIAL / READOUT
j0        = round(N/2);
j_sensors = (j0-99):(j0+100);  % 200 sensors centred on the stimulus

% Pre-calculate physical units for axes labeling
x_deg = linspace(-dov/2, dov/2, N);

n_type = numel(stim_types);
n_dur  = numel(dur_ms_list);

fprintf('\n==== RESULTS TIME-REVERSAL SYMMETRY BREAKING (1D) ====\n');
fprintf('Readout: 200 sensors\n');

ratio_result = zeros(n_type, n_dur);
V_result     = zeros(n_type, n_dur);
I_result     = zeros(n_type, n_dur);

figure(1)
set(gcf, 'Position', [50, 100, 320*n_dur, 260*n_type]);

%% RUN
for is = 1:n_type
    stype   = stim_types{is};
    v_dps   = v_list(is);
    per_deg = period_list(is);

    fprintf('\n--- %s (v=%.2f deg/s, period=%.2f deg) ---\n', stype, v_dps, per_deg);
    fprintf('  %-10s %-14s %-14s %-14s %-14s %-10s\n', ...
        'dur[ms]', 'R_fwd', 'R_rev', 'V', 'I', '|V|/|I|');

    for id = 1:n_dur
        dur_ms   = dur_ms_list(id);
        n_frames = round(dur_ms / 1000 * framerate);
        t_sec    = (0:n_frames-1) / framerate;

        S1D_fwd = make_TimeReversal_stimulus_1D(stype, N, n_frames, ...
            framerate, dov, phase0_rad, contrast, meanL, v_dps, per_deg);
        S1D_rev = flipud(S1D_fwd);   % time reversal, not direction flip

        O_fwd = stINRF3_compute_1D(S1D_fwd, dov, model, framerate, N, j_sensors);
        O_rev = stINRF3_compute_1D(S1D_rev, dov, model, framerate, N, j_sensors);

        R_fwd = stINRF3_readout_gaussian_1D(O_fwd, j_sensors, sigma_c);
        R_rev = stINRF3_readout_gaussian_1D(O_rev, j_sensors, sigma_c);
        V = (R_fwd - R_rev) / 2;
        I = (R_fwd + R_rev) / 2;
        ratio = abs(V) / max(abs(I), eps);

        V_result(is,id)     = V;
        I_result(is,id)     = I;
        ratio_result(is,id) = ratio;

        fprintf('  %-10d %-14.4e %-14.4e %-14.4e %-14.4e %-10.2f\n', ...
            dur_ms, R_fwd, R_rev, V, I, ratio);

        % Display and save the forward stimulus once per (type, duration)
        subplot(n_type, n_dur, (is-1)*n_dur + id)
        imagesc(x_deg, t_sec, S1D_fwd)
        colormap(gray)
        axis square
        xlabel('Space (degrees)')
        ylabel('Time (seconds)')
        title(sprintf('%s, %d ms', stype, dur_ms), 'Interpreter', 'none')
        pause(0.01)

        if make_videos
        video_filename = sprintf('TimeReversal_%s_dur%dms_fwd.avi', stype, dur_ms);
        v = VideoWriter(video_filename, 'Motion JPEG AVI');
        v.FrameRate = 60;
        v.Quality   = 100;
        open(v);

        height_2d = 400;
        S_normalized = S1D_fwd - min(S1D_fwd(:));
        S_uint8 = uint8(255 * (S_normalized / max(S_normalized(:))));
        for tt = 1:n_frames
            frame_2d = repmat(S_uint8(tt, :), height_2d, 1);
            writeVideo(v, frame_2d);
        end
        close(v);
        fprintf('Video saved: %s\n', video_filename);
        end
    end
end

fprintf('\n');

%% PLOT
fig = figure(2);
set(fig, 'Color', 'w', 'Position', [100 100 720 520]);
colors = [0.35 0.55 0.85; 0.85 0.30 0.45];
hold on; grid on; box on;
for is = 1:n_type
    plot(dur_ms_list, ratio_result(is,:), 'o-', ...
        'Color', colors(is,:), 'LineWidth', 1.8, ...
        'MarkerSize', 8, 'MarkerFaceColor', colors(is,:), ...
        'DisplayName', sprintf('%s (model)', stim_types{is}));
    yline(ratio_fly_l(is), '--', sprintf('%s fly = %.2f', stim_types{is}, ratio_fly_l(is)), ...
        'Color', colors(is,:), 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left', ...
        'Interpreter', 'none');
end
hold off;
xlabel('Duration (ms)', 'FontSize', 12);
ylabel('|V| / |I|', 'FontSize', 12);
title('Time-reversal symmetry breaking — model vs fly reference', ...
    'FontSize', 12, 'FontWeight', 'normal');
legend('Location', 'best', 'Interpreter', 'none');
set(gca, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out');

%% SAVE
outdir = fullfile(pwd, 'results_TimeReversal_1D');
if ~exist(outdir, 'dir'), mkdir(outdir); end

stamp = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(outdir, sprintf('TimeReversal_1D_%s.mat', stamp)), ...
    'ratio_result', 'V_result', 'I_result', 'stim_types', 'dur_ms_list', ...
    'ratio_fly_l', 'v_list', 'period_list', 'N', 'dov', 'framerate', ...
    'meanL', 'contrast', 'model');

saveas(fig, fullfile(outdir, sprintf('TimeReversal_1D_%s.png', stamp)));
fprintf('Saved results to: %s\n', outdir);
