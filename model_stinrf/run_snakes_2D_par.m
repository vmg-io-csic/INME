%% run_snakes_2D_par.m
%
% Parallel version of run_snakes_2D.m.
%
% Inputs
%   N              (int)    image side length [px]
%   dov            (double) field of view [dva]
%   framerate      (int)    temporal sampling [fps]
%   onset_ms       (int)    pattern onset time in ms
%   duration_ms    (int)    pattern visible duration in ms
%   post_offset_ms (int)    frames after offset in ms
%   mean_luminance (double) background luminance [a.u.]
%   L_illusion     (1x4 double) asymmetric luminance profile per quarter-period
%   L_symmetric    (1x4 double) symmetric control profile
%   orientations_deg (1xK double) pattern orientations to test in degrees
%   n_trials       (int)    repetitions per condition
%   rng_seed       (int)    random seed, [] for unseeded
%   half_region    (int)    half-side of the readout ROI in pixels
%   readout_mode   (char)   'single' | 'multi' (flat ROI) | 'gauss' (Gaussian c_s)
%   phase_mode     (char)   'fixed' or 'random'
%   vel_dps        (double) real-motion drift speed in degrees/second
%

clear; clc; close all;

cfg   = stINRF3_config('rotating_snakes');
model = cfg.model;

%% STIMULUS PARAMETERS
N              = 300;             % image size [px]
dov            = 300 * (1.5 / 128); % field of view [dva]
framerate      = 200;             % temporal sampling [fps]
onset_ms       = 150;             % pattern onset time [ms]
duration_ms    = 300;             % pattern visible duration [ms]
post_offset_ms = 200;             % frames after offset [ms]
mean_luminance = 50;              % background luminance [a.u.]

L_illusion     = [0,   5, 100, 50];
L_symmetric    = [0, 100,   0, 100];

orientations_deg = [0, 90];

%% TRIAL / READOUT
n_trials     = 5;
rng_seed     = 42;               % [] = no seed, or number
half_region  = 100;
readout_mode = 'gauss';          % 'single' | 'multi' | 'gauss'
phase_mode   = 'random';         % 'fixed' | 'random'
vel_dps      = 0.1;              % real_motion drift speed [deg/s]
                                 % Lower = smaller magnitude difference

%% DERIVED
total_ms  = onset_ms + duration_ms + post_offset_ms;
n_frames  = round(total_ms  / 1000 * framerate);
onset_f   = round(onset_ms  / 1000 * framerate);
offset_f  = round((onset_ms + duration_ms) / 1000 * framerate);
dpi       = N / dov;

if isfield(model, 'sf0_factor')
    sf0_cpd = model.sf0_factor / dov;
else
    sf0_cpd = model.sf0;
end
period_px = dpi / sf0_cpd;

if ~isempty(rng_seed), rng(rng_seed); end

%% CONDITIONS
cond_table = { ...
    'illusion',    L_illusion,  'moving'; ...
    'symmetric',   L_symmetric, 'zero';   ...
    'real_motion', L_illusion,  'moving'; ...
    'baseline',    [],          'zero';   ...
};

fprintf('\n==== RESULTS ROTATING SNAKES ====\n');
fprintf('Readout: %s | Phase: %s', readout_mode, phase_mode);
if strcmp(phase_mode, 'random'), fprintf(' (%d trials)', n_trials); end
fprintf('\n');

for ci = 1:size(cond_table,1)
    cond_name = cond_table{ci,1};
    L_table   = cond_table{ci,2};
    expected  = cond_table{ci,3};

    for oi = 1:numel(orientations_deg)
        ort = orientations_deg(oi);

        fprintf('\n--- %s | orientation: %d° ---\n', cond_name, ort);

        if strcmp(phase_mode, 'fixed')
            fprintf('  %-8s %-14s %-12s %-12s %-12s %-12s\n', ...
                'Trial', 'Phase', 'Theoric', 'Predicted', 'Error', 'Magnitude');

            TCB = make_stim(cond_name, N, n_frames, onset_f, offset_f, ...
                period_px, ort, L_table, L_illusion, mean_luminance, ...
                0, vel_dps, dpi, framerate);
            [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
            theta = do_readout(O_H, O_V, readout_mode, half_region);
            [~, mag] = do_readout_mag(O_H, O_V, readout_mode, half_region);

            if strcmp(expected, 'moving')
                err = stINRF3_angular_error(theta, ort);
                fprintf('  %-8d %-14s %-12s %-12s %-12s %-12s\n', 1, '0.0000°', ...
                    sprintf('%+.2f°', ort), sprintf('%+.2f°', theta), ...
                    sprintf('%+.2f°', err), sprintf('%.4f', mag));
            else
                fprintf('  %-8d %-14s %-12s %-12s %-12s %-12s\n', 1, '0.0000°', ...
                    'n/a', sprintf('%+.2f°', theta), 'n/a', sprintf('%.4f', mag));
            end

        else  % random
            fprintf('  %-8s %-14s %-12s %-12s %-12s %-12s\n', ...
                'Trial', 'Phase', 'Theoric', 'Predicted', 'Error', 'Magnitude');

            phases_px = rand(1, n_trials) * period_px;
            theta_all = zeros(1, n_trials);
            mag_all   = zeros(1, n_trials);

            % Compute all trials in parallel (each uses a pre-drawn phase)
            parfor tr = 1:n_trials
                TCB = make_stim(cond_name, N, n_frames, onset_f, offset_f, ...
                    period_px, ort, L_table, L_illusion, mean_luminance, ...
                    phases_px(tr), vel_dps, dpi, framerate);
                [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
                theta_all(tr) = do_readout(O_H, O_V, readout_mode, half_region);

                [~, mag_all(tr)] = do_readout_mag(O_H, O_V, readout_mode, half_region);
            end

            % Print results in order
            for tr = 1:n_trials
                if strcmp(expected, 'moving')
                    err_tr = stINRF3_angular_error(theta_all(tr), ort);
                    fprintf('  %-8d %-14s %-12s %-12s %-12s %-12s\n', tr, ...
                        sprintf('%.4f px', phases_px(tr)), ...
                        sprintf('%+.2f°', ort), sprintf('%+.2f°', theta_all(tr)), ...
                        sprintf('%+.2f°', err_tr), sprintf('%.4f', mag_all(tr)));
                else
                    fprintf('  %-8d %-14s %-12s %-12s %-12s %-12s\n', tr, ...
                        sprintf('%.4f px', phases_px(tr)), ...
                        'n/a', sprintf('%+.2f°', theta_all(tr)), ...
                        'n/a', sprintf('%.4f', mag_all(tr)));
                end
            end
            [m, s] = stINRF3_circular_mean(theta_all(:));
            fprintf('  Mean: %+.2f° (±%.2f°)  mag=%.4f\n', m, s, mean(mag_all));
        end
    end
end

fprintf('\n');


%% LOCAL FUNCTIONS
function TCB = make_stim(cond_name, N, n_frames, onset_f, offset_f, ...
        period_px, ort, L_table, L_illusion, mean_luminance, ...
        ph_px, vel_dps, dpi, framerate)
    switch cond_name
        case {'illusion','symmetric'}
            TCB = make_snakes_static(N, n_frames, onset_f, offset_f, ...
                period_px, ort, L_table, mean_luminance, ph_px);
        case 'real_motion'
            TCB = make_snakes_moving(N, n_frames, onset_f, offset_f, ...
                period_px, ort, L_illusion, mean_luminance, ...
                vel_dps, dpi, framerate, ph_px);
        case 'baseline'
            TCB = single(mean_luminance * ones(N, N, n_frames));
    end
end

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
function [d, mag] = do_readout_mag(O_H, O_V, mode, half_region)
% DO_READOUT_MAG  As do_readout, also returning the response magnitude.
    switch mode
        case 'single'
            [d, mag] = stINRF3_readout_single(O_H, O_V);
        case 'multi'
            [d, mag, ~] = stINRF3_readout_multisensor(O_H, O_V, half_region);
        case 'gauss'
            sigma_c = 2*half_region / sqrt(12);
            [d, mag, ~] = stINRF3_readout_gaussian(O_H, O_V, sigma_c, ...
                              'SupportHalf', half_region);
        otherwise
            error('do_readout_mag: unknown mode ''%s''.', mode);
    end
end

