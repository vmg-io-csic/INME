%% run_TimeReversal_trace_1D.m
%
% Response time courses for the time-reversal stimuli, averaged over spatial
% phase.
%
% Two traces are shown per condition: the raw instantaneous readout R(t),
% and the Tb-weighted contribution Tb(t)*R(t), whose sum is exactly the
% scalar reported by the sweep.
%
% NOTE ON PHASE AVERAGING. The spatial phase is randomised on every
% presentation. For the light-edge stimulus this changes very little here,
% because the reset time is set by mod(v*t, period) and does not depend on
% phase, and because the readout window spans an integer number of spatial
% periods. What it does buy is an error bar on every point.
%
% Inputs
%   ops (struct array) stype, Nm and Nw for each panel
%   dur_ms (double) stimulus duration [ms]
%   preview_ms (double) span of the space-time preview [ms]
%   n_phases (int) number of spatial phases to average
%   phase_mode (char) 'equispaced' (exact phase average) | 'random'
%   sigma_c (double) width of the Gaussian readout filter c_s [px]
%

clear; clc; close all;

cfg   = stINRF3_config('timereversal');
model = cfg.model;

ops(1).stype = 'LightEdges';        ops(1).Nm = 200;  ops(1).Nw = 300;
ops(2).stype = 'UnevenSquareWave';  ops(2).Nm = 750;  ops(2).Nw = 15;

N          = cfg.stim.N;
dov        = cfg.stim.dov;
framerate  = cfg.stim.framerate;
meanL      = cfg.stim.meanL;
contrast   = cfg.stim.contrast;

dur_ms     = 1000;
preview_ms = 250;
sigma_c    = 57.74;

n_phases   = 25;
phase_mode = 'equispaced';   % 'equispaced' | 'random'
rng_seed   = 42;

out_dir = stinrf_results_dir(mfilename);
if isempty(out_dir), out_dir = pwd; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

n_frames  = round(dur_ms / 1000 * framerate);
t_ms      = (0:n_frames-1)'/framerate*1000;
j0        = round(N/2);
j_sensors = (j0-99):(j0+100);
x_deg     = (0:N-1) * dov / N;
n_ops     = numel(ops);

switch phase_mode
    case 'equispaced'
        phases = 2*pi * (0:n_phases-1) / n_phases;
    case 'random'
        if ~isempty(rng_seed), rng(rng_seed); end
        phases = rand(1, n_phases) * 2*pi;
    otherwise
        error('phase_mode must be ''equispaced'' or ''random''.');
end

% Tb, the same beta window the engine applies internally
tn = (0:n_frames-1)' / n_frames;
Tb = tn.^model.lpa .* (1 - tn).^model.lpb;
Tb = Tb / (sum(Tb) + eps);

v_of   = struct('LightEdges', cfg.stim.v_LightEdges, ...
                'UnevenSquareWave', cfg.stim.v_UnevenSquareWave);
per_of = struct('LightEdges', cfg.stim.period_LightEdges, ...
                'UnevenSquareWave', cfg.stim.period_UnevenSquareWave);
fly_of = struct('LightEdges', cfg.stim.ratio_fly_LightEdges, ...
                'UnevenSquareWave', cfg.stim.ratio_fly_UnevenSquareWave);

R_fwd_t   = zeros(n_frames, n_ops);   R_rev_t   = zeros(n_frames, n_ops);
R_fwd_sem = zeros(n_frames, n_ops);   R_rev_sem = zeros(n_frames, n_ops);
R_fwd_per_phase    = zeros(n_phases, n_ops);
R_rev_per_phase    = zeros(n_phases, n_ops);
R_fwd_tb_per_phase = zeros(n_phases, n_ops);
R_rev_tb_per_phase = zeros(n_phases, n_ops);
stim_fwd  = cell(n_ops,1);  stim_rev = cell(n_ops,1);  rows = cell(n_ops,1);

fprintf('\n==== TIME-REVERSAL RESPONSE TRACES ====\n');
fprintf('%d ms at %d fps | readout Gaussian c_s, sigma_c = %.2f px\n', ...
        dur_ms, framerate, sigma_c);
fprintf('%d %s spatial phases\n', n_phases, phase_mode);

for k = 1:n_ops
    m = model;  m.Nm = ops(k).Nm;  m.Nw = ops(k).Nw;
    stype = ops(k).stype;
    v_k   = v_of.(stype);   per_k = per_of.(stype);

    Rf_ph = zeros(n_frames, n_phases);
    Rr_ph = zeros(n_frames, n_phases);

    parfor p = 1:n_phases
        S_fwd = make_TimeReversal_stimulus_1D(stype, N, n_frames, ...
            framerate, dov, phases(p), contrast, meanL, v_k, per_k);
        S_rev = flipud(S_fwd);

        [~, F_fwd] = stINRF3_compute_1D(S_fwd, dov, m, framerate, N, j_sensors);
        [~, F_rev] = stINRF3_compute_1D(S_rev, dov, m, framerate, N, j_sensors);

        rf = zeros(n_frames,1);  rr = zeros(n_frames,1);
        for tt = 1:n_frames
            rf(tt) = stINRF3_readout_gaussian_1D(F_fwd(tt,:), j_sensors, sigma_c);
            rr(tt) = stINRF3_readout_gaussian_1D(F_rev(tt,:), j_sensors, sigma_c);
        end
        Rf_ph(:,p) = rf;  Rr_ph(:,p) = rr;
    end

    R_fwd_per_phase(:,k)    = mean(Rf_ph, 1).';    % flat mean, one per phase
    R_rev_per_phase(:,k)    = mean(Rr_ph, 1).';
    R_fwd_tb_per_phase(:,k) = sum(Tb .* Rf_ph, 1).';
    R_rev_tb_per_phase(:,k) = sum(Tb .* Rr_ph, 1).';

    R_fwd_t(:,k)   = mean(Rf_ph, 2);
    R_rev_t(:,k)   = mean(Rr_ph, 2);
    R_fwd_sem(:,k) = std(Rf_ph, 0, 2) / sqrt(n_phases);
    R_rev_sem(:,k) = std(Rr_ph, 0, 2) / sqrt(n_phases);

    Rf_tb_ph = sum(Tb .* Rf_ph, 1);   Rr_tb_ph = sum(Tb .* Rr_ph, 1);
    Rf_fl_ph = mean(Rf_ph, 1);        Rr_fl_ph = mean(Rr_ph, 1);

    Rf_tb = mean(Rf_tb_ph);  Rr_tb = mean(Rr_tb_ph);
    Rf_fl = mean(Rf_fl_ph);  Rr_fl = mean(Rr_fl_ph);

    S0 = make_TimeReversal_stimulus_1D(stype, N, n_frames, framerate, ...
             dov, phases(1), contrast, meanL, v_k, per_k);
    stim_fwd{k} = S0;  stim_rev{k} = flipud(S0);

    rows{k} = struct('stype', stype, 'Nm', ops(k).Nm, 'Nw', ops(k).Nw, ...
        'R_fwd_tb', Rf_tb, 'R_rev_tb', Rr_tb, ...
        'V_tb', (Rf_tb-Rr_tb)/2, 'I_tb', (Rf_tb+Rr_tb)/2, ...
        'R_fwd_flat', Rf_fl, 'R_rev_flat', Rr_fl, ...
        'V_flat', (Rf_fl-Rr_fl)/2, 'I_flat', (Rf_fl+Rr_fl)/2, ...
        'sem_fwd_tb', std(Rf_tb_ph)/sqrt(n_phases), ...
        'sem_rev_tb', std(Rr_tb_ph)/sqrt(n_phases), ...
        'sem_fwd_flat', std(Rf_fl_ph)/sqrt(n_phases), ...
        'sem_rev_flat', std(Rr_fl_ph)/sqrt(n_phases), ...
        'fly_ratio', fly_of.(stype));

    fprintf('\n--- %s ---\n', disp_name(stype));
    fprintf('  %-14s %13s %13s %13s %13s %10s\n', ...
            'mean', 'R_fwd', 'R_rev', 'V', 'I', '|V|/|I|');
    fprintf('  %-14s %13.4e %13.4e %13.4e %13.4e %10.2f\n', 'Tb-weighted', ...
            Rf_tb, Rr_tb, rows{k}.V_tb, rows{k}.I_tb, ...
            abs(rows{k}.V_tb)/max(abs(rows{k}.I_tb), eps));
    fprintf('  %-14s %13.4e %13.4e\n', '   SEM', ...
            rows{k}.sem_fwd_tb, rows{k}.sem_rev_tb);
    fprintf('  %-14s %13.4e %13.4e %13.4e %13.4e %10.2f\n', 'flat', ...
            Rf_fl, Rr_fl, rows{k}.V_flat, rows{k}.I_flat, ...
            abs(rows{k}.V_flat)/max(abs(rows{k}.I_flat), eps));
    fprintf('  %-14s %13.4e %13.4e\n', '   SEM', ...
            rows{k}.sem_fwd_flat, rows{k}.sem_rev_flat);
end

V_t = (R_fwd_t - R_rev_t) / 2;
I_t = (R_fwd_t + R_rev_t) / 2;
R_fwd_tb_t = Tb .* R_fwd_t;
R_rev_tb_t = Tb .* R_rev_t;

stamp = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(out_dir, ['timereversal_trace_' stamp '.mat']), ...
     'ops','dur_ms','framerate','sigma_c','model','t_ms','x_deg','Tb', ...
     'n_phases','phase_mode','phases', ...
     'R_fwd_t','R_rev_t','R_fwd_sem','R_rev_sem', ...
     'R_fwd_per_phase','R_rev_per_phase', ...
     'R_fwd_tb_per_phase','R_rev_tb_per_phase', ...
     'R_fwd_tb_t','R_rev_tb_t','V_t','I_t','rows','stim_fwd','stim_rev');

col_fwd = [0.20 0.40 0.75];
col_rev = [0.80 0.25 0.20];
n_prev  = min(n_frames, round(preview_ms/1000*framerate));

figure('Color','w','Position',[60 60 900 460*n_ops]);
for k = 1:n_ops
    r1 = (2*k-2)*3;

    % --- stimulus, time forward
    subplot(2*n_ops, 3, [r1+1, r1+4]);
    imagesc(x_deg, t_ms(1:n_prev), stim_fwd{k}(1:n_prev,:));
    colormap(gray); set(gca,'YDir','reverse');
    xlabel('space (deg)'); ylabel('time (ms)');
    title('time forward','Color',col_fwd,'FontWeight','normal');
    set(gca,'XColor',col_fwd,'YColor',col_fwd,'LineWidth',1.5);

    % --- stimulus, time reversed
    subplot(2*n_ops, 3, [r1+2, r1+5]);
    imagesc(x_deg, t_ms(1:n_prev), stim_rev{k}(1:n_prev,:));
    colormap(gray); set(gca,'YDir','reverse');
    xlabel('space (deg)');
    title('time reversed','Color',col_rev,'FontWeight','normal');
    set(gca,'XColor',col_rev,'YColor',col_rev,'LineWidth',1.5);

    % --- response trace, time running downwards
    subplot(2*n_ops, 3, r1+3); hold on;
    yf = R_fwd_tb_t(:,k);  yr = R_rev_tb_t(:,k);
    sf = Tb .* R_fwd_sem(:,k);  sr = Tb .* R_rev_sem(:,k);
    fill([yf+sf; flipud(yf-sf)], [t_ms; flipud(t_ms)], col_fwd, ...
         'FaceAlpha', 0.25, 'EdgeColor', 'none');
    fill([yr+sr; flipud(yr-sr)], [t_ms; flipud(t_ms)], col_rev, ...
         'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(yf, t_ms, 'Color', col_fwd, 'LineWidth', 1.2);
    plot(yr, t_ms, 'Color', col_rev, 'LineWidth', 1.2);
    xline(0, 'k--', 'LineWidth', 0.5);
    set(gca, 'YDir', 'reverse');
    xlabel('R(t) (a.u.)'); ylabel('time (ms)');
    title(disp_name(ops(k).stype), 'FontWeight', 'bold');
    box on;

    % --- mean response
    subplot(2*n_ops, 3, r1+6); hold on;
    vals = [rows{k}.R_fwd_tb, rows{k}.R_rev_tb];
    sems = [rows{k}.sem_fwd_tb, rows{k}.sem_rev_tb];
    b = barh(1:2, vals, 0.5, 'FaceColor', 'flat');
    b.CData(1,:) = col_fwd;  b.CData(2,:) = col_rev;
    errorbar(vals, 1:2, sems, 'horizontal', 'k', 'LineStyle','none', 'LineWidth', 1);
    xline(0, 'k--', 'LineWidth', 0.5);
    set(gca, 'YTick', 1:2, 'YTickLabel', {'forward','reversed'}, 'YDir','reverse');
    xlabel('mean R (a.u.)');
    box on;
end

fprintf('\nSaved to %s\n\n', out_dir);

function s = disp_name(stype)
% DISP_NAME  Human-readable label for each stimulus.
    switch stype
        case 'LightEdges',       s = 'Light Edges';
        case 'UnevenSquareWave', s = 'Uneven Square Wave';
        otherwise,               s = stype;
    end
end
