%% check_readout_gaussian.m
%
% Compares the square ROI against a sweep of Gaussian readout widths on
% gratings and Type-II plaids, evaluating every window on the same O maps.
%
% Inputs
%   sigmas       (1xS double) sigma_c values to sweep [px]
%   support_half (int)        half-side of the square support [px]
%   n_phases     (int)        random-phase repetitions per stimulus
%

clear; clc;

sigmas       = [16.7 25 33.33 40 50 57.74 70];  % [px]
support_half = 100;         % keep the legacy 200x200 support
n_phases     = 5;
rng(42);

%% ---------------- shared stimulus / model parameters -------------------
N           = 300;
dov         = 300*(1.5/128);
framerate   = 200;
duration_ms = 2000;
f_cpd       = 2.0;
speed       = 8 / framerate;
contrast    = 1;  meanL = 100;
n_frames    = round(duration_ms/1000*framerate);
f_cpx       = f_cpd / (N/dov);

cfg   = stINRF3_config('gratings');
model = cfg.model;

fprintf('sigma_c in px | dpi = %.2f px/deg | carrier period = %.1f px\n', ...
        N/dov, 1/f_cpx);

%% ---------------- 1) gratings -----------------------------------------
for dir_deg = [0 30 45]
    phases = rand(1, n_phases) * 2*pi;
    O = cell(n_phases, 2);
    for k = 1:n_phases
        TCB = make_grating(N, n_frames, f_cpx, speed, dir_deg, phases(k), contrast, meanL);
        [O{k,1}, O{k,2}] = stINRF3_compute(TCB, dov, model, framerate, N);
    end
    report(sprintf('grating %d deg', dir_deg), O, sigmas, support_half);
end

%% ---------------- 2) plaid II -----------------------------------------
% configs: [speed1, speed2, theta1, theta2, ioc]
configs = [1 2 -45 -25  27.12
           1 2  50  35 -25.95
           1 2 -50 -30  22.12];
for c = 1:size(configs,1)
    th1 = configs(c,3);  th2 = configs(c,4);
    beat_cpd = 2*f_cpd*abs(sind((th1-th2)/2));
    fprintf('\n[plaid II cfg %d]  beat = %.3f cpd -> period %.1f px\n', ...
            c, beat_cpd, (N/dov)/beat_cpd);
    ph = rand(n_phases, 2) * 2*pi;
    O  = cell(n_phases, 2);
    for k = 1:n_phases
        TCB = make_plaidII(N, n_frames, f_cpx, configs(c,1)*speed, configs(c,2)*speed, ...
                           th1, th2, ph(k,1), ph(k,2), contrast, meanL);
        [O{k,1}, O{k,2}] = stINRF3_compute(TCB, dov, model, framerate, N);
    end
    report(sprintf('plaid II cfg %d (IOC %+.2f)', c, configs(c,5)), O, sigmas, support_half);
end

%% ---------------- helper ----------------------------------------------
function report(name, O, sigmas, support_half)
    n = size(O,1);
    sq = zeros(n,1);
    for k = 1:n
        sq(k) = stINRF3_readout_multisensor(O{k,1}, O{k,2}, support_half);
    end
    fprintf('\n=== %s ===\n', name);
    fprintf('  %-22s %10s %10s %12s\n', 'readout', 'mean', 'spread', 'vs square');
    fprintf('  %-22s %10.4f %10.4f %12s\n', ...
            sprintf('square %dx%d', 2*support_half, 2*support_half), ...
            circmean(sq), max(sq)-min(sq), '-');
    for s = sigmas
        g = zeros(n,1);
        for k = 1:n
            g(k) = stINRF3_readout_gaussian(O{k,1}, O{k,2}, s, ...
                       'SupportHalf', support_half);
        end
        fprintf('  %-22s %10.4f %10.4f %12.4f\n', ...
                sprintf('gauss sigma_c=%.2f', s), ...
                circmean(g), max(g)-min(g), circmean(g)-circmean(sq));
    end
end

function m = circmean(a)
    a = deg2rad(a(:));
    m = rad2deg(atan2(mean(sin(a)), mean(cos(a))));
end
