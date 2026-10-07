%% run_Interaction_1D.m
%
% Cross-scale interaction between a low and a high spatial frequency
% component, swept over
% the five composition conditions and the width of the Gaussian temporal
% envelope.
%
% Conditions:
%   1 S1  — 1 cpd alone, moving
%   2 S3  — 3 cpd alone, moving
%   3 CS  — 1 cpd moving  + 3 cpd static
%   4 CSr — 3 cpd moving  + 1 cpd static   (direction reversal expected
%           at short sigma_t)
%   5 CM  — both moving
%
% Sign convention: a positive response means the model reports rightward
% motion, negative means leftward. A trial counts as correct when the sign
% of the response matches the sign of `direction`.
%
% NOTE ON THE FRAME RATE. 120 fps is deliberate, not a placeholder: the
% effect does not appear at 1000 fps with this parameter set. Since the
% model integrates over Nw = 40 ms, the frame rate sets how many samples
% fall inside that window (5 frames at 120 fps, 40 at 1000), so it is part
% of the operating point rather than a sampling detail.
%
% Inputs
%   N            (int)    number of spatial pixels
%   dov          (double) field of view [dva]
%   framerate    (int)    frame rate [fps]
%   duration_ms  (int)    total stimulus duration [ms]
%   speed_dps    (double) drift speed [degrees/second]
%   meanL        (double) mean luminance [a.u.]
%   sfm_Low_cpd  (double) low component spatial frequency [cpd]
%   sfm_High_cpd (double) high component spatial frequency [cpd]
%   sigma_t_list (1xS double) Gaussian envelope widths to sweep [ms]
%   n_trials     (int)    random-phase repetitions per cell
%   rng_seed     (int)    random seed, 0 for unseeded
%   sigma_c      (double) width of the Gaussian readout filter c_s [px]
%

clear; clc; close all;

cfg   = stINRF3_config('DH');
model = cfg.model;

%% STIMULUS PARAMETERS
N            = 400;          % spatial samples [px]
dov          = 4.0;          % field of view [dva]
framerate    = 1000;         % temporal sampling [fps]
duration_ms  = 500;          % stimulus duration [ms]
motion_mode  = 'tf';         % 'tf'    = both components drift at tf_hz Hz
                             %           (the Derrington & Henning convention)
                             % 'speed' = both drift at speed_dps deg/s
tf_hz        = 12;           % temporal frequency [Hz], used when mode is 'tf'
speed_dps    = 8;            % drift speed [deg/s], used when mode is 'speed'
meanL        = 100;          % mean luminance [a.u.]
sfm_Low_cpd  = 1;            % low component spatial frequency [cpd]
sfm_High_cpd = 3;            % high component spatial frequency [cpd]
direction    = -1;           % -1 leftward (as characterised), +1 rightward

% Full range, spanning the reversal, with a fine grid across the crossing.
% The transition is abrupt -- 0 per cent at 87.5 ms and above 90 at 90 ms --
% so the crossing cannot be found by interpolating between widely spaced
% widths: the straight line that interpolation assumes is not there. The
% widths between 86 and 90 ms are simulated instead, and the figure takes
% whichever of them measured closest to chance.
sigma_t_list = [10, 12.5, 25, 37.5, 50, 62.5, 75, 85, 86, 87, 88, 88.5, ...
                89, 89.2, 89.5, 90, 95, 100, 125, 150, 200];

% Per-component drift speeds implied by the motion mode
switch motion_mode
    case 'tf'
        v_Low  = tf_hz / sfm_Low_cpd;
        v_High = tf_hz / sfm_High_cpd;
    case 'speed'
        v_Low  = speed_dps;
        v_High = speed_dps;
    otherwise
        error('run_Interaction_1D: motion_mode must be ''tf'' or ''speed''.');
end

cond_all_names  = {'S1', 'S3', 'CS', 'CSr', 'CM'};
cond_all_labels = {'1cpd moving', '3cpd moving', '1mov + 3stat', ...
                   '3mov + 1stat', 'both moving'};

% Which conditions to run. Use 4 alone to hunt the CSr threshold quickly.
cond_subset = 1:5;
cond_names  = cond_all_names(cond_subset);
cond_labels = cond_all_labels(cond_subset);

%% TRIAL
n_trials = 25;               % random-phase repetitions
% Seeded. Unseeded, the crossing moved by about a millisecond from run to
% run and any width quoted for it went stale the next time the experiment
% was run. Note that with 25 trials the share correct steps by four points,
% so it can reach 48 or 52 per cent but never exactly 50.
rng_seed = 1;                % 0 = unseeded (different each run)

%% READOUT
sigma_c = 57.74;             % [px]

%% DERIVED
n_frames  = round(duration_ms / 1000 * framerate);
n_cond    = numel(cond_names);
n_sigma   = numel(sigma_t_list);
j0        = round(N/2);
j_sensors = (j0-99):(j0+100);   % 200 sensors centred on the stimulus
if rng_seed > 0, rng(rng_seed); end

% Pre-draw every phase so the parfor is reproducible and order-independent
ph_high = rand(n_cond, n_sigma, n_trials) * 2*pi;
ph_low  = rand(n_cond, n_sigma, n_trials) * 2*pi;

% Flatten condition x sigma_t x trial into a single parfor loop
n_jobs = n_cond * n_sigma * n_trials;
[ic_g, is_g, it_g] = ndgrid(1:n_cond, 1:n_sigma, 1:n_trials);
ic_vec = ic_g(:);  is_vec = is_g(:);  it_vec = it_g(:);

fprintf('\n==== RESULTS CROSS-SCALE INTERACTION (1D) ====\n');
if direction > 0, dir_str = 'rightward'; else, dir_str = 'leftward'; end
fprintf('%d conditions x %d envelope widths x %d trials | %d ms | %d fps | motion %s\n', ...
        n_cond, n_sigma, n_trials, duration_ms, framerate, dir_str);
fprintf('Motion: %s | low %.2f cpd at %.2f deg/s (%.1f Hz) | high %.2f cpd at %.2f deg/s (%.1f Hz)\n', ...
        motion_mode, sfm_Low_cpd, v_Low, sfm_Low_cpd*v_Low, ...
        sfm_High_cpd, v_High, sfm_High_cpd*v_High);
fprintf('Readout: Gaussian c_s, sigma_c = %.2f px\n', sigma_c);

%% RUN
val_flat = zeros(1, n_jobs);

parfor jj = 1:n_jobs
    ic = ic_vec(jj);  is = is_vec(jj);  it = it_vec(jj);

    mI = make_Interaction_stimulus_1D(N, n_frames, framerate, dov, ...
        cond_subset(ic), sigma_t_list(is), ph_high(ic,is,it), ph_low(ic,is,it), ...
        v_High, v_Low, sfm_High_cpd, sfm_Low_cpd, meanL, direction);

    O_sensors = stINRF3_compute_1D(mI, dov, model, framerate, N, j_sensors);
    val_flat(jj) = stINRF3_readout_gaussian_1D(O_sensors, j_sensors, sigma_c);
end

val = reshape(val_flat, n_cond, n_sigma, n_trials);

%% SAVE (before printing, so a formatting error cannot lose the run)
out_dir = stinrf_results_dir(mfilename);
if isempty(out_dir), out_dir = pwd; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
save(fullfile(out_dir, ['interaction_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), ...
     'cond_names', 'cond_labels', 'cond_subset', 'sigma_t_list', 'duration_ms', ...
     'framerate', 'direction', 'motion_mode', 'tf_hz', 'speed_dps', ...
     'v_Low', 'v_High', 'sigma_c', 'model', 'val');

%% RESULTS
pct_correct = 100 * mean(val * direction > 0, 3);
mean_resp   = mean(val, 3);
sem_resp    = std(val, 0, 3) / sqrt(n_trials);

fprintf('\n--- Sign matches stimulus direction (%%) ---\n');
fprintf('  %-6s %-16s', 'Cond', 'Composition');
for is = 1:n_sigma, fprintf('%12s', sprintf('%g ms', sigma_t_list(is))); end
fprintf('\n  %s\n', repmat('-', 1, 23 + 12*n_sigma));
for ic = 1:n_cond
    fprintf('  %-6s %-16s', cond_names{ic}, cond_labels{ic});
    fprintf('%11.0f%%', pct_correct(ic,:));
    fprintf('\n');
end

fprintf('\n--- Mean response (a.u.) ---\n');
fprintf('  %-6s %-16s', 'Cond', 'Composition');
for is = 1:n_sigma, fprintf('%12s', sprintf('%g ms', sigma_t_list(is))); end
fprintf('\n  %s\n', repmat('-', 1, 23 + 12*n_sigma));
for ic = 1:n_cond
    fprintf('  %-6s %-16s', cond_names{ic}, cond_labels{ic});
    fprintf('%12.4f', mean_resp(ic,:));
    fprintf('\n');
end

fprintf('\n--- SEM of the response (a.u.) ---\n');
fprintf('  %-6s %-16s', 'Cond', 'Composition');
for is = 1:n_sigma, fprintf('%12s', sprintf('%g ms', sigma_t_list(is))); end
fprintf('\n  %s\n', repmat('-', 1, 23 + 12*n_sigma));
for ic = 1:n_cond
    fprintf('  %-6s %-16s', cond_names{ic}, cond_labels{ic});
    fprintf('%12.4f', sem_resp(ic,:));
    fprintf('\n');
end

%% THRESHOLD: sigma_t at which performance crosses chance (50%)
% Linear interpolation in log(sigma_t) on the first upward crossing of 50%.
fprintf('\n--- Chance-level crossing (50%% correct) ---\n');
for ic = 1:n_cond
    p = pct_correct(ic,:);
    k = find(p(1:end-1) < 50 & p(2:end) >= 50, 1, 'last');
    if isempty(k)
        if all(p >= 50)
            fprintf('  %-6s always above chance\n', cond_names{ic});
        elseif all(p < 50)
            fprintf('  %-6s always below chance\n', cond_names{ic});
        else
            fprintf('  %-6s no clean crossing (non-monotonic)\n', cond_names{ic});
        end
        continue;
    end
    l1 = log(sigma_t_list(k));  l2 = log(sigma_t_list(k+1));
    thr = exp(l1 + (50 - p(k)) * (l2 - l1) / (p(k+1) - p(k)));
    fprintf('  %-6s sigma_t = %6.1f ms   (between %g and %g ms)\n', ...
            cond_names{ic}, thr, sigma_t_list(k), sigma_t_list(k+1));
end

% Same crossing read from the mean response changing sign
fprintf('\n--- Sign reversal of the mean response ---\n');
for ic = 1:n_cond
    m = mean_resp(ic,:) * direction;   % >0 means correct direction
    k = find(m(1:end-1) < 0 & m(2:end) >= 0, 1, 'last');
    if isempty(k)
        fprintf('  %-6s no sign reversal across the sweep\n', cond_names{ic});
        continue;
    end
    l1 = log(sigma_t_list(k));  l2 = log(sigma_t_list(k+1));
    thr = exp(l1 + (0 - m(k)) * (l2 - l1) / (m(k+1) - m(k)));
    fprintf('  %-6s sigma_t = %6.1f ms   (between %g and %g ms)\n', ...
            cond_names{ic}, thr, sigma_t_list(k), sigma_t_list(k+1));
end
fprintf('\n');

%% FIGURE
figure('Color','w','Position',[100 100 640 460]); hold on;
markers = {'o','s','d','^','v'};
for ic = 1:n_cond
    errorbar(sigma_t_list, mean_resp(ic,:), sem_resp(ic,:), ...
        ['-' markers{ic}], 'LineWidth', 1.5, 'MarkerSize', 7, ...
        'MarkerFaceColor', 'w');
end
yline(0, 'k--', 'LineWidth', 1);
set(gca, 'XScale', 'log');
xlabel('Temporal envelope \sigma_t (ms)');
ylabel('Mean response (a.u.)');
legend(cond_names, 'Location', 'best');
title('Cross-scale interaction: response vs envelope width');
grid on; box on;


export_Interaction;
