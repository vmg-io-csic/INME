%% run_plaidI_2D.m
%
% Runs the Type-I plaid experiment using the 2D model.
%
% Inputs
%   N              (int)    spatial resolution [px]
%   dov            (double) field of view in degrees of visual angle
%   framerate      (int)    frame rate [fps]
%   duration_ms    (int)    stimulus duration in ms
%   f_spatial_cpd  (double) spatial frequency per component in cycles/degree
%   speed_Hz       (double) temporal frequency in Hz
%   contrast       (double) Michelson contrast per component
%   meanL          (double) mean luminance
%   angle_offset   (double) component offset from rigid direction in degrees
%   rigid_dirs_deg (1xD double) rigid directions to sweep in degrees
%   n_trials       (int)    repetitions per direction
%   rng_seed       (int)    random seed, [] for unseeded
%   half_region    (int)    half-side of the readout ROI in pixels
%   readout_mode   (char)   'single' | 'multi' (flat ROI) | 'gauss' (Gaussian c_s)
%   phase_mode     (char)   'fixed' or 'random'
%

clear; clc; close all;

cfg   = stINRF3_config('plaidI');
model = cfg.model;

%% STIMULUS PARAMETERS
N              = 300;              % image size [px]
dov            = 300*(1.5/128);    % field of view [dva]
framerate      = 200;              % temporal sampling [fps]
duration_ms    = 2000;             % stimulus duration [ms]
f_spatial_cpd  = 2.0;             % spatial frequency per component [cpd]
speed_Hz       = 8;               % temporal frequency [Hz]
speed          = speed_Hz / framerate;   % [cy/frame]
contrast       = 1;               % Michelson contrast per component
meanL          = 100;             % mean luminance [a.u.]
angle_offset   = 45;              % component offset from rigid direction [deg]
rigid_dirs_deg = 0:15:90;         % rigid directions to sweep [deg]

%% TRIAL / READOUT
n_trials     = 5;                 % random-phase repetitions
rng_seed     = 42;                % random seed ([] = unseeded)
half_region  = 100;               % ROI half-side [px]
readout_mode = 'gauss';           % 'single' | 'multi' | 'gauss'
phase_mode   = 'random';          % 'fixed' | 'random'

%% DERIVED
n_dir    = numel(rigid_dirs_deg);
n_frames = round(duration_ms / 1000 * framerate);
f_cpx    = f_spatial_cpd / (N / dov);

if ~isempty(rng_seed), rng(rng_seed); end

fprintf('\n==== RESULTS PLAID I ====\n');
fprintf('Readout: %s | Phase: %s', readout_mode, phase_mode);
if strcmp(phase_mode, 'random'), fprintf(' (%d trials)', n_trials); end
fprintf('\n');

target = rigid_dirs_deg(:);

if strcmp(phase_mode, 'fixed')
    for id = 1:n_dir
        t1 = rigid_dirs_deg(id)+angle_offset;  t2 = rigid_dirs_deg(id)-angle_offset;
        fprintf('\n--- Direction: %d° ---\n', rigid_dirs_deg(id));
        fprintf('  %-8s %-24s %-12s %-12s %-12s\n', 'Trial', 'Phases', 'Theoric', 'Predicted', 'Error');

        TCB = make_plaidI(N, n_frames, f_cpx, speed, t1, t2, 0, 0, contrast, meanL);
        [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
        pred_d = do_readout(O_H, O_V, readout_mode, half_region);
        err_d  = stINRF3_angular_error(pred_d, rigid_dirs_deg(id));

        fprintf('  %-8d %-24s %-12s %-12s %-12s\n', 1, ...
            'phi1=0.0000, phi2=0.0000', ...
            sprintf('%+.2f°', rigid_dirs_deg(id)), ...
            sprintf('%+.2f°', pred_d), ...
            sprintf('%+.2f°', err_d));
    end

else
    phases1 = rand(n_dir, n_trials) * 2*pi;
    phases2 = rand(n_dir, n_trials) * 2*pi;

    % Flatten direction x trial into a single parfor loop, so the pool can
    % use up to n_dir*n_trials workers instead of n_trials.
    n_jobs = n_dir * n_trials;
    id_vec = repelem(1:n_dir, n_trials);
    tr_vec = repmat(1:n_trials, 1, n_dir);
    d_flat = zeros(1, n_jobs);

    parfor jj = 1:n_jobs
        id = id_vec(jj);  tr = tr_vec(jj);
        t1 = rigid_dirs_deg(id) + angle_offset;
        t2 = rigid_dirs_deg(id) - angle_offset;
        TCB = make_plaidI(N, n_frames, f_cpx, speed, t1, t2, ...
            phases1(id,tr), phases2(id,tr), contrast, meanL);
        [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
        d_flat(jj) = do_readout(O_H, O_V, readout_mode, half_region);
    end

    d_mat = reshape(d_flat, n_trials, n_dir);

    for id = 1:n_dir
        fprintf('\n--- Direction: %d° ---\n', rigid_dirs_deg(id));
        fprintf('  %-8s %-24s %-12s %-12s %-12s\n', 'Trial', 'Phases', 'Theoric', 'Predicted', 'Error');

        d_all = d_mat(:, id);
        for tr = 1:n_trials
            err_tr = stINRF3_angular_error(d_all(tr), rigid_dirs_deg(id));
            ph_str = sprintf('phi1=%.4f, phi2=%.4f', rad2deg(phases1(id,tr)), rad2deg(phases2(id,tr)));
            fprintf('  %-8d %-24s %-12s %-12s %-12s\n', tr, ph_str, ...
                sprintf('%+.2f°', rigid_dirs_deg(id)), ...
                sprintf('%+.2f°', d_all(tr)), ...
                sprintf('%+.2f°', err_tr));
        end
        [m, s] = stINRF3_circular_mean(d_all);
        fprintf('  Mean: %+.2f° (±%.2f°)\n', m, s);
    end
end

fprintf('\n');

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
