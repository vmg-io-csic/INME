%% run_LightEdges.m
%
% Data for the summary figure of Nm=200, Nw=300, Tb(a=1,b=1), dur=1000ms
% (V/I=2.18, I/SEM=36.4 -- MUCH more significant than the Nm=75,Nw=150
% Tb(9,1) combo, found in check_lightedges_phases.m/log). Reproduces the
% EXACT same 20 phases (rng_seed=7, direct rand(n_phases,1) draw, no
% causality-trick reweighting -- dur=1000ms is simulated directly) so
% the bar chart matches that log's own numbers exactly.
%
% n_phases is 20, the number the sweep was read at. rand(n,1) after a fixed
% seed is a prefix of one stream, so 20 phases here are the first 20 of any
% longer draw: the earlier 40-phase export averaged those same 20 plus 20
% more and reported V/I=2.15 rather than 2.18.

clear; clc;

cfg = stINRF3_config('timereversal');
model = cfg.model;
model.lpa = 1; model.lpb = 1;
model.Nm = 200; model.Nw = 300;

N          = cfg.stim.N;
dov        = cfg.stim.dov;
framerate  = cfg.stim.framerate;
meanL      = cfg.stim.meanL;
contrast   = cfg.stim.contrast;
v_dps      = cfg.stim.v_LightEdges;
period_deg = cfg.stim.period_LightEdges;

j0        = round(N/2);
j_sensors = (j0-99):(j0+100);
sigma_c   = 57.74;

dur_ms   = 1000;
nf       = round(dur_ms/1000*framerate);
n_phases = 20;
rng_seed = 7;

if ~isempty(rng_seed), rng(rng_seed); end
phases = rand(n_phases,1) * 2*pi;

out_dir = fullfile(fileparts(mfilename('fullpath')), 'results_best_lightedges2');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% Stimulus space-time matrices (first of the 20 phases, for display)
S_stim  = make_TimeReversal_stimulus_1D('LightEdges', N, nf, framerate, dov, phases(1), contrast, meanL, v_dps, period_deg);
Sr_stim = flipud(S_stim);
writematrix(S_stim,  fullfile(out_dir, 'stim_fwd.csv'));
writematrix(Sr_stim, fullfile(out_dir, 'stim_rev.csv'));

%% R(t) trace, same representative phase, raw (instantaneous)
[~, F_fwd] = stINRF3_compute_1D(S_stim, dov, model, framerate, N, j_sensors);
[~, F_rev] = stINRF3_compute_1D(Sr_stim, dov, model, framerate, N, j_sensors);

r_fwd_raw = zeros(nf,1); r_rev_raw = zeros(nf,1);
for tt = 1:nf
    r_fwd_raw(tt) = stINRF3_readout_gaussian_1D(F_fwd(tt,:), j_sensors, sigma_c);
    r_rev_raw(tt) = stINRF3_readout_gaussian_1D(F_rev(tt,:), j_sensors, sigma_c);
end
t = (0:nf-1)' * 1000/framerate;

T = table(t, r_fwd_raw, r_rev_raw, 'VariableNames', {'t_ms','R_fwd_raw','R_rev_raw'});
writetable(T, fullfile(out_dir, 'trace.csv'));

%% Bar chart: the SAME 40 phases as check_lightedges_phases.m, Tb-weighted mean
Rf_bar = zeros(n_phases,1); Rr_bar = zeros(n_phases,1);
for p = 1:n_phases
    Sp_fwd = make_TimeReversal_stimulus_1D('LightEdges', N, nf, framerate, dov, phases(p), contrast, meanL, v_dps, period_deg);
    Sp_rev = flipud(Sp_fwd);
    O_fwd = stINRF3_compute_1D(Sp_fwd, dov, model, framerate, N, j_sensors);
    O_rev = stINRF3_compute_1D(Sp_rev, dov, model, framerate, N, j_sensors);
    Rf_bar(p) = stINRF3_readout_gaussian_1D(O_fwd, j_sensors, sigma_c);
    Rr_bar(p) = stINRF3_readout_gaussian_1D(O_rev, j_sensors, sigma_c);
    if mod(p,5)==0, fprintf('  phase %d/%d done\n', p, n_phases); end
end
Tbar = table((1:n_phases)', Rf_bar, Rr_bar, 'VariableNames', {'phase_idx','R_fwd','R_rev'});
writetable(Tbar, fullfile(out_dir, 'bar.csv'));

Imean = (mean(Rf_bar)+mean(Rr_bar))/2; Isem = std((Rf_bar+Rr_bar)/2)/sqrt(n_phases);
fprintf('mean R_fwd=%.4f  mean R_rev=%.4f\n', mean(Rf_bar), mean(Rr_bar));
fprintf('I=%.4f (SEM=%.4f, |I|/SEM=%.2f)  V=%.4f  V/I=%.4f\n', Imean, Isem, abs(Imean)/Isem, ...
    (mean(Rf_bar)-mean(Rr_bar))/2, (mean(Rf_bar)-mean(Rr_bar))/(mean(Rf_bar)+mean(Rr_bar)));
fprintf('All data written to %s\n', out_dir);

export_LightEdges;
