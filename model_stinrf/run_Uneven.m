%% run_Uneven.m
%
% Data for the summary figure of Nm=40, Nw=10, Tb(a=1,b=1), dur=1000ms
% (V/I=-12.27, the closest ratio match to the -12.74 target of anything
% tested, though NOT statistically significant: |I|/SEM=1.10).
%
% Reproduces the EXACT same 20 phases that sweep_uneven_fine2.m used for
% this combo (same rng(42) call, same phases_all(n_combos,n_phases)
% construction, same loop-order index) so the bar chart matches the
% sweep's own numbers exactly, instead of drawing a fresh phase set.
%
% n_phases has to be 20 and not merely at least 20. The draw here is a
% matrix, rand(n_combos, n_phases), so its width is part of the stream:
% asking for 25 or 40 phases does not add phases to row 21, it deals the
% whole matrix differently and row 21 becomes a different set. An earlier
% 40-phase export did exactly that and reported V/I = -45 for this combo.

clear; clc;

cfg = stINRF3_config('timereversal');
model = cfg.model;
model.lpa = 1; model.lpb = 1;
model.Nm = 40; model.Nw = 10;

N          = cfg.stim.N;
dov        = cfg.stim.dov;
framerate  = cfg.stim.framerate;
meanL      = cfg.stim.meanL;
contrast   = cfg.stim.contrast;
v_dps      = cfg.stim.v_UnevenSquareWave;
period_deg = cfg.stim.period_UnevenSquareWave;

j0        = round(N/2);
j_sensors = (j0-99):(j0+100);
sigma_c   = 57.74;

dur_ms   = 1000;
n_frames = round(dur_ms/1000*framerate);
n_phases = 20;
rng_seed = 42;

% Reconstruct sweep_uneven_fine2.m's exact phase matrix and index
Nm_list = [20 25 30 35 40];
Nw_list = [10 12 15 18 20];
n_combos = numel(Nm_list) * numel(Nw_list);
rng(rng_seed);
phases_all = rand(n_combos, n_phases) * 2*pi;
idx = (find(Nm_list==40)-1)*numel(Nw_list) + find(Nw_list==10);  % = 21
phases = phases_all(idx, :);
fprintf('combo index = %d\n', idx);

out_dir = fullfile(fileparts(mfilename('fullpath')), 'results_best_uneven3');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% Stimulus space-time matrices (first of the 20 phases, for display)
S_stim  = make_TimeReversal_stimulus_1D('UnevenSquareWave', N, n_frames, framerate, dov, phases(1), contrast, meanL, v_dps, period_deg);
Sr_stim = flipud(S_stim);
writematrix(S_stim,  fullfile(out_dir, 'stim_fwd.csv'));
writematrix(Sr_stim, fullfile(out_dir, 'stim_rev.csv'));

%% R(t) trace, same representative phase, raw (instantaneous)
[~, F_fwd] = stINRF3_compute_1D(S_stim, dov, model, framerate, N, j_sensors);
[~, F_rev] = stINRF3_compute_1D(Sr_stim, dov, model, framerate, N, j_sensors);

r_fwd_raw = zeros(n_frames,1); r_rev_raw = zeros(n_frames,1);
for tt = 1:n_frames
    r_fwd_raw(tt) = stINRF3_readout_gaussian_1D(F_fwd(tt,:), j_sensors, sigma_c);
    r_rev_raw(tt) = stINRF3_readout_gaussian_1D(F_rev(tt,:), j_sensors, sigma_c);
end
t = (0:n_frames-1)' * 1000/framerate;

T = table(t, r_fwd_raw, r_rev_raw, 'VariableNames', {'t_ms','R_fwd_raw','R_rev_raw'});
writetable(T, fullfile(out_dir, 'trace.csv'));

%% Bar chart: the SAME 20 phases as the sweep, Tb-weighted mean
Rf_bar = zeros(n_phases,1); Rr_bar = zeros(n_phases,1);
for p = 1:n_phases
    Sp_fwd = make_TimeReversal_stimulus_1D('UnevenSquareWave', N, n_frames, framerate, dov, phases(p), contrast, meanL, v_dps, period_deg);
    Sp_rev = flipud(Sp_fwd);
    O_fwd = stINRF3_compute_1D(Sp_fwd, dov, model, framerate, N, j_sensors);
    O_rev = stINRF3_compute_1D(Sp_rev, dov, model, framerate, N, j_sensors);
    Rf_bar(p) = stINRF3_readout_gaussian_1D(O_fwd, j_sensors, sigma_c);
    Rr_bar(p) = stINRF3_readout_gaussian_1D(O_rev, j_sensors, sigma_c);
end
Tbar = table((1:n_phases)', Rf_bar, Rr_bar, 'VariableNames', {'phase_idx','R_fwd','R_rev'});
writetable(Tbar, fullfile(out_dir, 'bar.csv'));

Imean = (mean(Rf_bar)+mean(Rr_bar))/2; Isem = std((Rf_bar+Rr_bar)/2)/sqrt(n_phases);
fprintf('mean R_fwd=%.4f  mean R_rev=%.4f\n', mean(Rf_bar), mean(Rr_bar));
fprintf('I=%.4f (SEM=%.4f, |I|/SEM=%.2f)  V=%.4f  V/I=%.4f\n', Imean, Isem, abs(Imean)/Isem, ...
    (mean(Rf_bar)-mean(Rr_bar))/2, (mean(Rf_bar)-mean(Rr_bar))/(mean(Rf_bar)+mean(Rr_bar)));
fprintf('All data written to %s\n', out_dir);

export_Uneven;
