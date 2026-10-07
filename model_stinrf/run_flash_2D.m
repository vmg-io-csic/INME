%% run_flash_2D.m
%
% Runs the static flash (ON-OFF disc transient) experiment.
%
% Inputs
%   N               (int)    image size [px], square
%   dov             (double) field of view [dva]
%   framerate       (int)    frame rate [fps]
%   duration_ms     (int)    total duration [ms]
%   flash_onset_ms  (int)    disc onset time in ms
%   flash_dur_ms    (int)    disc visible duration in ms
%   mean_luminance  (double) background luminance [a.u.]
%   disc_amplitude  (double) disc luminance increment [a.u.]
%   disc_radius_deg (double) disc radius in degrees
%   ring_width_px   (int)    annulus half-width for ring readout in pixels
%   n_arrows        (int)    number of direction arrows around disc edge
%
clear; clc; close all;

cfg   = stINRF3_config('flash');
model = cfg.model;

%% STIMULUS PARAMETERS
N               = 128;   % image size [px]
dov             = 1.5;   % field of view [dva]
framerate       = 1000;  % temporal sampling [fps]
duration_ms     = 750;   % total stimulus duration [ms]
flash_onset_ms  = 150;   % disc onset time [ms]
flash_dur_ms    = 400;   % disc visible duration [ms]
mean_luminance  = 100;   % background luminance [a.u.]
disc_amplitude  = 50;    % disc luminance increment [a.u.]
disc_radius_deg = 0.3;   % disc radius [dva]

%% READOUT PARAMETERS
ring_width_px = 5;       % annulus half-width for ring readout [px]
n_arrows      = 20;      % direction arrows around disc edge

%% DERIVED
n_frames        = round(duration_ms / 1000 * framerate);
dpi             = N / dov;
disc_r_px       = round(disc_radius_deg * dpi);
flash_offset_ms = flash_onset_ms + flash_dur_ms;

%% MASKS
disc_mask = make_disc_mask(N, dov, disc_radius_deg);

cx = (N+1)/2;  cy = (N+1)/2;
[XX, YY]  = meshgrid(1:N, 1:N);
dist_map  = sqrt((XX - cx).^2 + (YY - cy).^2);

mask_ring     = dist_map >= (disc_r_px - ring_width_px) & ...
                dist_map <= (disc_r_px + ring_width_px);
mask_interior = dist_map <  (disc_r_px - ring_width_px);
mask_exterior = dist_map >  (disc_r_px + ring_width_px);

%% ARROW POSITIONS
ang_arr  = linspace(0, 2*pi*(1 - 1/n_arrows), n_arrows);
arr_x_px = round(cx + disc_r_px * cos(ang_arr));
arr_y_px = round(cy + disc_r_px * sin(ang_arr));
arr_x_px = max(1, min(N, arr_x_px));
arr_y_px = max(1, min(N, arr_y_px));

%% GENERATE STIMULUS
TCB_sq = make_flash(N, n_frames, framerate, disc_mask, ...
    flash_onset_ms, flash_dur_ms, mean_luminance, disc_amplitude, 'square');

%% RUN MODEL WITH TEMPORAL OUTPUT
[O_H, O_V, O_H_t, O_V_t] = stINRF3_compute(TCB_sq, dov, model, framerate, N);

%% KEY FRAMES
t_onset_show  = flash_onset_ms  + model.Nm;
t_mid         = flash_onset_ms  + round(flash_dur_ms / 2);
t_offset_show = flash_offset_ms + model.Nm;
t_end         = min(flash_offset_ms + 4*model.Nm, duration_ms - 1);

key_times_ms  = [0, t_onset_show, t_mid, t_offset_show, t_end];
key_labels    = {'Baseline', ...
    sprintf('Onset (t=%dms)', t_onset_show), ...
    sprintf('Mid-flash (t=%dms)', t_mid), ...
    sprintf('Offset (t=%dms)', t_offset_show), ...
    sprintf('End (t=%dms)', t_end)};
key_frames    = min(key_times_ms + 1, n_frames);
n_keys        = numel(key_frames);

%% PER-PIXEL MAGNITUDE
mag_t = sqrt(O_H_t.^2 + O_V_t.^2);

%% RING SCALAR TIME SERIES
ring_scalar = zeros(1, n_frames);
for i = 1:n_frames
    fm = mag_t(:,:,i);
    ring_scalar(i) = mean(fm(mask_ring));
end

%% ============================================================
%%  PRINT RESULTS
%% ============================================================
fprintf('\n==== RESULTS STATIC FLASH ====\n\n');

% --- Per-frame ring magnitude ---
fprintf('--- Per-frame magnitude by zone ---\n');
fprintf('  %-22s %12s %12s %12s\n', 'Frame', 'Ring', 'Interior', 'Exterior');
fprintf('  %s\n', repmat('-', 1, 60));

for k = 1:n_keys
    fi = key_frames(k);
    fm = mag_t(:,:,fi);
    fprintf('  %-22s %12.4f %12.4f %12.4f\n', key_labels{k}, ...
        mean(fm(mask_ring)), mean(fm(mask_interior)), mean(fm(mask_exterior)));
end

% --- Direction at disc edge ---
for k = 2:n_keys  % skip baseline
    fi = key_frames(k);
    oh_frame = O_H_t(:,:,fi);
    ov_frame = O_V_t(:,:,fi);
    fm = mag_t(:,:,fi);
    ring_mag_k = mean(fm(mask_ring));

    fprintf('\n--- Direction at disc edge (%d positions, %s) ---\n', n_arrows, key_labels{k});
    fprintf('  ');
    for a = 1:n_arrows
        ah = oh_frame(arr_y_px(a), arr_x_px(a));
        av = ov_frame(arr_y_px(a), arr_x_px(a));
        am = sqrt(ah^2 + av^2);
        if am > 0.01 * ring_mag_k
            fprintf('%+4.0f° ', atan2d(av, ah));
        else
            fprintf('  ·  ');
        end
    end
    fprintf('\n');
end

% --- Vectorial cancellation (Tb-integrated) ---
mag_total = sqrt(O_H.^2 + O_V.^2);
scalar_ring = mean(mag_total(mask_ring));
vec_oh = mean(O_H(mask_ring));
vec_ov = mean(O_V(mask_ring));
vec_mag = sqrt(vec_oh^2 + vec_ov^2);

fprintf('\n--- Vectorial cancellation (Tb-integrated, ring) ---\n');
fprintf('  Scalar mag (mean per-pixel):   %.4f\n', scalar_ring);
fprintf('  Vectorial mag (mean O_H,O_V):  %.4e\n', vec_mag);
fprintf('  Cancellation ratio:            %.0f×\n', scalar_ring / (vec_mag + eps));

%% SAVE .mat for plotting
save('results_flash_2D.mat', ...
    'model', 'N', 'dov', 'framerate', 'duration_ms', ...
    'flash_onset_ms', 'flash_dur_ms', 'mean_luminance', 'disc_amplitude', ...
    'disc_radius_deg', 'disc_r_px', 'ring_width_px', 'n_arrows', ...
    'O_H', 'O_V', 'O_H_t', 'O_V_t', 'mag_t', ...
    'ring_scalar', 'key_times_ms', 'key_labels', 'key_frames', ...
    'mask_ring', 'mask_interior', 'mask_exterior', ...
    'arr_x_px', 'arr_y_px', 'ang_arr', 'cx', 'cy', ...
    'TCB_sq', 'disc_mask');
fprintf('\nSaved: results_flash_2D.mat (for plot_flash.m)\n');
