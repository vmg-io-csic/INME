%% run_gratings_2D.m
%
% Runs the single sinusoidal grating experiment using the 2D model.
%
% Inputs
%   N              (int)    image size [px] (square)
%   dov            (double) field of view in degrees of visual angle
%   framerate      (int)    sampling rate [fps]
%   duration_ms    (int)    stimulus duration in ms
%   f_spatial_cpd  (double) spatial frequency in cycles/degree
%   speed_Hz       (double) temporal frequency in Hz
%   contrast       (double) Michelson contrast, typically 1
%   meanL          (double) mean luminance [a.u.]
%   directions_deg (1xD double) grating directions to sweep in degrees
%   n_trials       (int)    repetitions per direction (random phase mode)
%   rng_seed       (int)    random seed, [] for unseeded
%   half_region    (int)    half-side of the readout ROI in pixels
%   readout_mode   (char)   'single' | 'multi' (flat ROI) | 'gauss' (Gaussian c_s)
%   phase_mode     (char)   'fixed' or 'random'
%

clear; clc; close all;

cfg   = stINRF3_config('gratings');
model = cfg.model;

%% STIMULUS PARAMETERS
N              = 300;
dov            = 300*(1.5/128);
framerate      = 200;
duration_ms    = 2000;
f_spatial_cpd  = 2.0;
speed_Hz       = 8;
speed          = speed_Hz / framerate;
contrast       = 1;
meanL          = 100;
directions_deg = 0:15:90;

%% TRIAL / READOUT
n_trials     = 5;
rng_seed     = 42;
half_region  = 100;
readout_mode = 'gauss';
phase_mode   = 'random';

%% DERIVED
n_dir    = numel(directions_deg);
n_frames = round(duration_ms / 1000 * framerate);
f_cpx    = f_spatial_cpd / (N / dov);

if ~isempty(rng_seed), rng(rng_seed); end

fprintf('\n==== RESULTS GRATINGS ====\n');
fprintf('Readout: %s | Phase: %s', readout_mode, phase_mode);
if strcmp(phase_mode, 'random'), fprintf(' (%d trials)', n_trials); end
fprintf('\n');

target = directions_deg(:);

if strcmp(phase_mode, 'fixed')
    for id = 1:n_dir
        fprintf('\n--- Direction: %d° ---\n', directions_deg(id));
        fprintf('  %-8s %-12s %-12s %-12s %-12s\n', 'Trial', 'Phase', 'Theoric', 'Predicted', 'Error');

        TCB = make_grating(N, n_frames, f_cpx, speed, directions_deg(id), 0, contrast, meanL);
        [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
        pred_d = do_readout(O_H, O_V, readout_mode, half_region);
        err_d  = stINRF3_angular_error(pred_d, directions_deg(id));

        fprintf('  %-8d %-12s %-12s %-12s %-12s\n', 1, '0.0000°', ...
            sprintf('%+.2f°', directions_deg(id)), ...
            sprintf('%+.2f°', pred_d), ...
            sprintf('%+.2f°', err_d));
    end

else
    phases = rand(n_dir, n_trials) * 2*pi;
    for id = 1:n_dir
        fprintf('\n--- Direction: %d° ---\n', directions_deg(id));
        fprintf('  %-8s %-12s %-12s %-12s %-12s\n', 'Trial', 'Phase', 'Theoric', 'Predicted', 'Error');

        d_all = zeros(n_trials, 1);
        for tr = 1:n_trials
            TCB = make_grating(N, n_frames, f_cpx, speed, directions_deg(id), phases(id,tr), contrast, meanL);
            [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
            d_all(tr) = do_readout(O_H, O_V, readout_mode, half_region);
            err_tr = stINRF3_angular_error(d_all(tr), directions_deg(id));

            fprintf('  %-8d %-12s %-12s %-12s %-12s\n', tr, ...
                sprintf('%.4f°', rad2deg(phases(id,tr))), ...
                sprintf('%+.2f°', directions_deg(id)), ...
                sprintf('%+.2f°', d_all(tr)), ...
                sprintf('%+.2f°', err_tr));
        end
        [m, s] = stINRF3_circular_mean(d_all);
        fprintf('  Mean: %+.2f° (±%.2f°)\n', m, s);
    end
end

fprintf('\n');

function d = do_readout(O_H, O_V, mode, half_region)
% Spatial readout.
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
