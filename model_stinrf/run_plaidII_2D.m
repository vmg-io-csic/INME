%% run_plaidII_2D.m
%
% Type-II plaids: the IOC and VA predictions differ, and human observers
% follow the IOC. Sweeps the IOC direction from -30 to +30 deg at the
% characterised stimulus operating point (2.5 cpd, components near 7 and
% 14 Hz).
%
% Inputs
%   N, dov, framerate, duration_ms, contrast, meanL  stimulus geometry
%   f_spatial_cpd (double) spatial frequency per component [cpd]
%   configs (Cx5) [speed1, speed2, theta1_deg, theta2_deg, ioc_dir_deg]
%   use_random_phases (logical) false = single run at phase 0. Full-field
%                     plaids are phase invariant under the Gaussian readout,
%                     so one run is normally enough.
%   n_trials (int) random-phase repetitions when the above is true
%   rng_seed (int) random seed, [] for unseeded
%   half_region (int) half-side of the readout ROI [px]
%   readout_mode (char) 'single' | 'multi' (flat ROI) | 'gauss' (Gaussian c_s)
%
clear; clc; close all;

cfg   = stINRF3_config('plaidII');
model = cfg.model;

%% STIMULUS PARAMETERS
N             = 300;              % image size [px]
dov           = 300*(1.5/128);    % field of view [dva]
framerate     = 200;              % temporal sampling [fps]
duration_ms   = 2000;             % stimulus duration [ms]
f_spatial_cpd = 2.5;              % spatial frequency per component [cpd]
contrast      = 1;                % Michelson contrast per component
meanL         = 100;              % mean luminance [a.u.]

%% CONFIGURATIONS
% Speed unit at f_spatial_cpd = 2.5: 1.0 corresponds to 11.72 Hz.
configs = [ ...
    0.6034, 1.2068,  -45,  -25,   27.12;  ... %  1   7.07 / 14.14 Hz
    0.6034, 1.2068,   50,   35,  -25.95;  ... %  2   7.07 / 14.14 Hz
    0.6034, 1.2068,  -50,  -30,   22.12;  ... %  3   7.07 / 14.14 Hz
    0.6046, 1.2045,   42,   22,  -30.00;  ... %  4   7.08 / 14.11 Hz
    0.6046, 1.2045,   52,   32,  -20.00;  ... %  5   7.08 / 14.11 Hz
    0.6661, 1.0931,   54,   34,  -10.00;  ... %  6   7.81 / 12.81 Hz
    0.7185, 1.0134,   54,   34,    0.00;  ... %  7   8.42 / 11.88 Hz
    0.6661, 1.0931,  -54,  -34,   10.00;  ... %  8   7.81 / 12.81 Hz
    0.6046, 1.2045,  -52,  -32,   20.00;  ... %  9   7.08 / 14.11 Hz
    0.6046, 1.2045,  -42,  -22,   30.00];    % 10   7.08 / 14.11 Hz

%% PHASE / READOUT
use_random_phases = true;
n_trials          = 5;
rng_seed          = 42;
half_region       = 100;
readout_mode      = 'gauss';

%% DERIVED
n_cfg    = size(configs, 1);
n_frames = round(duration_ms / 1000 * framerate);
f_cpx    = f_spatial_cpd / (N / dov);
tf_scale = f_spatial_cpd / 0.21333 / framerate;   % relative speed -> cy/frame

if use_random_phases, n_runs = n_trials; else, n_runs = 1; end
if ~isempty(rng_seed), rng(rng_seed); end

%% TARGETS: IOC and VA recomputed from the table itself
ioc_t = zeros(n_cfg,1);
va_t  = zeros(n_cfg,1);
t2ok  = false(n_cfg,1);
for c = 1:n_cfg
    [ioc_t(c), va_t(c), ~, t2ok(c)] = plaidII_ioc(configs(c,1), configs(c,2), ...
                                                  configs(c,3), configs(c,4));
end
if any(abs(stINRF3_angular_error(ioc_t, configs(:,5))) > 0.05)
    warning('run_plaidII_2D:iocMismatch', ...
        'Column 5 disagrees with the IOC direction recomputed from the table.');
end
if ~all(t2ok)
    warning('run_plaidII_2D:notTypeII', ...
        '%d config(s) are not Type-II.', sum(~t2ok));
end

fprintf('\n==== PLAID II ====\n');
if use_random_phases
    fprintf('Configs %d | %d random-phase runs each | readout %s\n', ...
        n_cfg, n_runs, readout_mode);
else
    fprintf('Configs %d | single run at phase 0 | readout %s\n', ...
        n_cfg, readout_mode);
end

%% RUN
phases1 = zeros(n_cfg, n_runs);
phases2 = zeros(n_cfg, n_runs);
if use_random_phases
    phases1 = rand(n_cfg, n_runs) * 2*pi;
    phases2 = rand(n_cfg, n_runs) * 2*pi;
end

pred_runs = zeros(n_cfg, n_runs);
pred      = zeros(n_cfg, 1);
sdev      = zeros(n_cfg, 1);

fprintf('\n--- Results ---\n');
fprintf('  %-6s %-8s %-8s %-8s %-8s %-9s %-9s %-9s %-10s\n', ...
    'Cfg', 'speed1', 'speed2', 'theta1', 'theta2', 'IOC_deg', 'VA_deg', ...
    'pred_deg', 'err_vs_IOC');

for c = 1:n_cfg
    s1  = configs(c,1) * tf_scale;
    s2  = configs(c,2) * tf_scale;
    th1 = configs(c,3);
    th2 = configs(c,4);

    d_all = zeros(n_runs, 1);
    parfor tr = 1:n_runs
        TCB = make_plaidII(N, n_frames, f_cpx, s1, s2, th1, th2, ...
            phases1(c,tr), phases2(c,tr), contrast, meanL);
        [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
        d_all(tr)  = do_readout(O_H, O_V, readout_mode, half_region);
    end
    pred_runs(c,:) = d_all(:).';

    if n_runs > 1
        [pred(c), sdev(c)] = stINRF3_circular_mean(d_all);
        for tr = 1:n_runs
            fprintf('  %-6s %-8.2f %-8.2f %-8.2f %-8.2f %-9.2f %-9.2f %-9.2f %-10.2f\n', ...
                sprintf('%d.%d', c, tr), configs(c,1), configs(c,2), ...
                configs(c,3), configs(c,4), ioc_t(c), va_t(c), d_all(tr), ...
                stINRF3_angular_error(d_all(tr), ioc_t(c)));
        end
        fprintf('  %-6d %-8s %-8s %-8s %-8s %-9s %-9s %-9.2f %-10.2f  (SD %.2f)\n\n', ...
            c, '', '', '', '', '', '', pred(c), ...
            stINRF3_angular_error(pred(c), ioc_t(c)), sdev(c));
    else
        pred(c) = d_all(1);
        fprintf('  %-6d %-8.2f %-8.2f %-8.2f %-8.2f %-9.2f %-9.2f %-9.2f %-10.2f\n', ...
            c, configs(c,1), configs(c,2), configs(c,3), configs(c,4), ...
            ioc_t(c), va_t(c), pred(c), ...
            stINRF3_angular_error(pred(c), ioc_t(c)));
    end
end

err_ioc = stINRF3_angular_error(pred, ioc_t);
err_va  = stINRF3_angular_error(pred, va_t);

%% SUMMARY
fprintf('\n--- Summary ---\n');
fprintf('  MAE vs IOC %.2f deg | MAE vs VA %.2f deg\n', ...
    mean(abs(err_ioc)), mean(abs(err_va)));
[worst, wc] = max(abs(err_ioc));
fprintf('  Worst config %d, error %.2f deg\n', wc, worst);
if n_runs > 1
    fprintf('  Mean circular SD across configs %.2f deg\n', mean(sdev));
end

%% FIGURE
figure('Color','w','Position',[100 100 560 460]); hold on;
lim = [-40 40];
plot(lim, lim, 'k--', 'LineWidth', 1.2);
plot(ioc_t, va_t, 'x', 'MarkerSize', 10, 'LineWidth', 1.2, ...
    'Color', [0.65 0.65 0.65]);
if n_runs > 1
    errorbar(ioc_t, pred, sdev, 'o', 'MarkerSize', 9, 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'w', 'CapSize', 4);
else
    plot(ioc_t, pred, 'o', 'MarkerSize', 9, 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'w');
end
xlim(lim); ylim([-90 90]);
xlabel('IOC direction (deg)'); ylabel('Predicted direction (deg)');
legend({'IOC prediction','VA prediction','stINRF3'}, 'Location','northwest');
title(sprintf('Plaid II: MAE %.2f deg vs IOC, %.2f deg vs VA', ...
    mean(abs(err_ioc)), mean(abs(err_va))));
grid on; box on; axis square;

%% LOCAL FUNCTIONS
function d = do_readout(O_H, O_V, mode, half_region)
% DO_READOUT  Spatial readout selector.
    switch mode
        case 'single'
            [d,~] = stINRF3_readout_single(O_H, O_V);
        case 'multi'
            [d,~,~] = stINRF3_readout_multisensor(O_H, O_V, half_region);
        case 'gauss'
            sigma_c = 2*half_region / sqrt(12);   % = 57.74 px for half_region = 100
            [d,~,~] = stINRF3_readout_gaussian(O_H, O_V, sigma_c, ...
                          'SupportHalf', half_region);
        otherwise
            error('do_readout: unknown mode ''%s''.', mode);
    end
end
