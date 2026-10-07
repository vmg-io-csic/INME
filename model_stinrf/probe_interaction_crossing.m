%% probe_interaction_crossing.m
%
% Finds the envelope width at which the CSr condition sits at chance.
%
% run_Interaction_1D reports the crossing by interpolating in log(sigma_t)
% between the two swept widths either side of 50 per cent, so the number it
% prints was never simulated. This script simulates the neighbourhood
% instead, on a fine grid, and reports what each width actually measured.
%
% Only CSr is run: it is the one condition that reverses, and the only one
% the figure uses. Fifty trials rather than twenty-five because with
% twenty-five the possible answers step by four points and land on 48 or 52
% but never on 50. The seed is fixed, so the answer is the same next time;
% the sweep in run_Interaction_1D leaves it at 0 and its crossing estimate
% moves by about a millisecond from run to run.
%
% Writes nothing. Reads its answer off the printed table.

clear; clc;

cfg   = stINRF3_config('DH');
model = cfg.model;

N            = 400;
dov          = 4.0;
framerate    = 1000;
duration_ms  = 500;
tf_hz        = 12;
meanL        = 100;
sfm_Low_cpd  = 1;
sfm_High_cpd = 3;
direction    = -1;

cond      = 4;                    % CSr: 3 cpd moving on a 1 cpd static
sigma_t_list = 85:0.5:90;
n_trials  = 50;
rng_seed  = 1;

v_Low  = tf_hz / sfm_Low_cpd;
v_High = tf_hz / sfm_High_cpd;

n_frames  = round(duration_ms / 1000 * framerate);
n_sigma   = numel(sigma_t_list);
j0        = round(N/2);
j_sensors = (j0-99):(j0+100);
sigma_c   = 57.74;

rng(rng_seed);
ph_high = rand(n_sigma, n_trials) * 2*pi;
ph_low  = rand(n_sigma, n_trials) * 2*pi;

n_jobs = n_sigma * n_trials;
[is_g, it_g] = ndgrid(1:n_sigma, 1:n_trials);
is_vec = is_g(:);  it_vec = it_g(:);

% The saved Processes profile caps at 8 workers and the machine has 32
% cores free, so the cap is raised on this cluster object only. Not saved,
% so the profile on disk is left as it was.
p = gcp('nocreate');
if isempty(p)
    c = parcluster('Processes');
    c.NumWorkers = 16;
    parpool(c, 16);
end

fprintf('\n==== CSr crossing probe ====\n');
fprintf('%d widths x %d trials | %d ms | %d fps | motion leftward | seed %d\n', ...
        n_sigma, n_trials, duration_ms, framerate, rng_seed);

val_flat = zeros(1, n_jobs);
parfor jj = 1:n_jobs
    is = is_vec(jj);  it = it_vec(jj);
    mI = make_Interaction_stimulus_1D(N, n_frames, framerate, dov, ...
        cond, sigma_t_list(is), ph_high(is,it), ph_low(is,it), ...
        v_High, v_Low, sfm_High_cpd, sfm_Low_cpd, meanL, direction);
    O = stINRF3_compute_1D(mI, dov, model, framerate, N, j_sensors);
    val_flat(jj) = stINRF3_readout_gaussian_1D(O, j_sensors, sigma_c);
end
val = reshape(val_flat, n_sigma, n_trials);

pct  = 100 * mean(val * direction > 0, 2);
mu   = mean(val, 2);
psem = 100 * sqrt((pct/100) .* (1 - pct/100) / n_trials);

fprintf('\n  sigma_t     correct        SEM    mean response\n');
fprintf('  ------------------------------------------------\n');
for is = 1:n_sigma
    fprintf('  %6.1f ms  %6.1f%%  +/- %4.1f%%    %+10.6f\n', ...
            sigma_t_list(is), pct(is), psem(is), mu(is));
end
[~, ib] = min(abs(pct - 50));
fprintf('\n  closest to chance: %.1f ms at %.1f%% correct\n', ...
        sigma_t_list(ib), pct(ib));
