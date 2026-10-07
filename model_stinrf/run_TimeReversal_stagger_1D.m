%% run_TimeReversal_stagger_1D.m
%
% Tests whether the sign oscillation in the light-edge response is caused
% by the synchronised global reset of the stimulus.
%
% All tiles resetting at once makes the whole field step from light to dark
% once per cycle, which drives the biphasic filter Tm and produces a large
% transient. Staggering the resets across space spreads that step out, so
% the field-average luminance stays nearly constant. If the hypothesis is
% right, the ratio |mean| / (peak-to-peak) of the response should rise
% sharply with n_tiles while the mean itself stays roughly the same.
%
% Also reports the symmetry signature of each stimulus, since staggering
% changes the unit cell: only variants that keep chi-Gamma-Theta (and
% nothing else) remain comparable.
%
% Inputs
%   n_tiles_list (1xK int) tiles per spatial period to sweep
%   Nm, Nw (double) temporal parameters of the model [ms]
%   dur_ms (double) stimulus duration [ms]
%   sigma_c (double) Gaussian readout width [px]
%

clear; clc; close all;

cfg   = stINRF3_config('timereversal');
model = cfg.model;

%% PARAMETERS
n_tiles_list = [1 2 3 6 12];
Nm      = 200;               % [ms]
Nw      = 300;               % [ms]
dur_ms  = 1000;              % [ms]
sigma_c = 57.74;             % [px]

N          = cfg.stim.N;
dov        = cfg.stim.dov;
framerate  = cfg.stim.framerate;
meanL      = cfg.stim.meanL;
contrast   = cfg.stim.contrast;
phase0_rad = cfg.stim.phase0_rad;
v_dps      = cfg.stim.v_LightEdges;
period_deg = cfg.stim.period_LightEdges;

model.Nm = Nm;  model.Nw = Nw;

%% DERIVED
n_frames  = round(dur_ms/1000*framerate);
t_ms      = (0:n_frames-1)'/framerate*1000;
j0        = round(N/2);
j_sensors = (j0-99):(j0+100);
nK        = numel(n_tiles_list);

fprintf('\n==== LIGHT EDGES: STAGGERED RESETS ====\n');
fprintf('Nm = %g ms, Nw = %g ms | %d ms at %d fps | sigma_c = %.2f px\n', ...
        Nm, Nw, dur_ms, framerate, sigma_c);
fprintf('reset frequency = %.2f Hz\n\n', v_dps/period_deg);

R_all   = zeros(n_frames, nK);
summary = zeros(nK, 5);

for k = 1:nK
    nt = n_tiles_list(k);
    [S, lum] = make_LightEdges_staggered_1D(N, n_frames, framerate, dov, ...
                   phase0_rad, contrast, meanL, v_dps, period_deg, nt);

    [~, F] = stINRF3_compute_1D(S, dov, model, framerate, N, j_sensors);
    R = zeros(n_frames,1);
    for tt = 1:n_frames
        R(tt) = stINRF3_readout_gaussian_1D(F(tt,:), j_sensors, sigma_c);
    end
    R_all(:,k) = R;

    lum_exc = max(lum) - min(lum);
    lum_jump = max(abs(diff(lum)));
    pk2pk   = max(R) - min(R);
    ratio   = 100*abs(mean(R))/max(pk2pk, eps);
    summary(k,:) = [nt, lum_exc, lum_jump, mean(R), ratio];

    fprintf('  n_tiles %2d | luminance excursion %.3f, max jump %.3f | mean R %+8.4f | peak-to-peak %8.4f | ratio %5.1f%%\n', ...
            nt, lum_exc, lum_jump, mean(R), pk2pk, ratio);
end

%% SYMMETRY SIGNATURE OF EACH VARIANT
% Built on exactly one spatial period and one temporal cycle so that the
% cyclic shifts can express the phase offsets.
fprintf('\n--- Symmetry signature (0 = symmetric) ---\n');
M    = 120;
fr_s = round(M / (period_deg/v_dps));       % one cycle in M frames
for k = 1:nK
    nt = n_tiles_list(k);
    Sp = make_LightEdges_staggered_1D(M, M, fr_s, period_deg, 0, 1, 100, ...
                                      v_dps, period_deg, nt);
    [d, names] = check_stimulus_symmetries(Sp);
    fprintf('  n_tiles %2d :', nt);
    for i = 1:numel(d)
        fprintf('  %s %.3f', names{i}, d(i));
    end
    if d(7) < 1e-6 && all(d(1:6) > 1e-6)
        fprintf('   <- chi-Gamma-Theta only');
    end
    fprintf('\n');
end

%% FIGURE
figure('Color','w','Position',[80 80 980 420]);
subplot(1,2,1); hold on;
for k = 1:nK
    plot(t_ms, R_all(:,k), 'LineWidth', 1.1);
end
yline(0,'k-','LineWidth',0.5);
xlabel('time (ms)'); ylabel('R(t) (a.u.)');
title(sprintf('Light edges, N_m=%g N_w=%g', Nm, Nw), 'Interpreter','tex');
legend(compose('%d tiles', n_tiles_list), 'Location','best');
grid on; box on;

subplot(1,2,2);
yyaxis left;  plot(summary(:,1), summary(:,2), '-o', 'LineWidth', 1.5);
ylabel('luminance excursion');
yyaxis right; plot(summary(:,1), summary(:,5), '-s', 'LineWidth', 1.5);
ylabel('|mean| / peak-to-peak (%)');
xlabel('tiles per period'); grid on; box on;
title('Global step vs response regularity');

fprintf('\n');
