%% run_CentralDrift_2D.m
%
% Runs the "Central Drift" illusion through stINRF3. The pattern is a ring of
% identical elements, each carrying a continuous asymmetric luminance gradient
% that runs from dark to light; it is built by make_CentralDrift_stimulus_2D.m
% from the per-element template in central_drift_template.m. Unlike Rotating
% Snakes, this
% illusion is designed to work under direct central fixation, not just
% peripherally -- hence "central" drift.
%
% Two conditions (baseline omitted -- a uniform grey canvas reliably
% gives ~0 in every prior 2D experiment here, not worth re-verifying):
%   forward   — the pattern as-is. Looks CW. Expected: positive.
%   invertido — average of the fliplr- and flipud-reflected image's
%               tangential-projection FIELDS (before smoothing), which
%               cancels the cardinal bias each single reflection
%               introduces in the model's H/V-only Gabor channels.
%               Looks CCW. Expected: negative.
%
% Read-out: global tangential projection (make_tangential_field.m +
% smooth_and_score_proj.m), sigma profile centred on the wheel (=canvas
% centre here, since this is a centred-only experiment), averaged over
% the full disc, EXACTLY as in run_RotatingSnakesWheel_2D.m -- no
% sign correction applied here. proj_raw is used as returned by
% make_tangential_field.m (positive = CW per that function's docstring,
% same formula, same convention, no per-stimulus adjustment).
%
% Verified against María Berraquero's direct visual check of the
% generated forward/invertido videos: 'forward' looks CW and 'invertido'
% looks CCW -- the SAME pairing as Rotating Snakes -- consistent with
% proj_tan's raw sign being positive for 'forward' in both stimuli.
% No sign adjustment needed anywhere in this script.
%
% IMPORTANT -- n_segments override: the template has 12 elements, but 12 is a
% multiple of 4, and this metric is EXACTLY
% zero for any n-fold pattern with n a multiple of 4 when the model only
% has H/V (0/90 deg) Gabor channels (the same degeneracy diagnosed and
% fixed for the Rotating Snakes wheel, n_units 24->23; verified here:
% forward with n_segments=12 gives proj_tan ~5e-11, not a real null).
% This script therefore generates the stimulus with n_segments=13 (one
% extra pill, barely visible) to break the degeneracy -- verified to
% give a real, non-degenerate signal (+1.63e-3 for 'forward').
%
% Model parameters: sf0_factor = N/100 = 3 (same 100 px Gabor period
% convention as the other 2D experiments); sigmaMu = 0.12 deg;
% sigmaW = 0.22 deg; Nw = 75 ms (all validated optima carried over from
% the Rotating Snakes wheel work -- not independently re-swept here).
%
% For each condition, the stimulus is displayed and saved as a video,
% and a final bar plot of the tangential-projection scalar is saved to
% results_CentralDrift_2D/.
%
% Inputs

%   N              (int)    canvas size [px], square
%   n_segments     (int)    pill count used for model runs (13, see above)
%   dov            (double) field of view [dva]
%   framerate      (int)    temporal sampling [fps]
%   onset_ms, duration_ms, post_off_ms (int) temporal paradigm [ms]
%   meanL          (double) background luminance [a.u.]
%   sigma_profile  (char)   smoothing profile, see smooth_and_score_proj.m
%
clear; clc; close all;

%% STIMULUS PARAMETERS
N           = 300;
n_segments  = 13;       % breaks the 4-fold exact-zero degeneracy; reference has 12
dov         = N * (1.5 / 128);
framerate   = 200;
onset_ms    = 150;
duration_ms = 150;
post_off_ms = 200;
meanL       = 50;
contrast    = 1.0;

sigma_profile = 'D';

%% MODEL PARAMETERS
cfg   = stINRF3_config('rotating_snakes');
model = cfg.model;
model.sf0    = [];
model = rmfield(model, 'sf0');
model.sf0_factor = round(N / 100);
model.sigmaMu    = 0.12;
model.sigmaW     = 0.22;
model.Nw         = 75;

%% READOUT
cx_wheel = (N+1)/2;
cy_wheel = (N+1)/2;
r_inner_ro = 0;

%% DERIVED
total_ms = onset_ms + duration_ms + post_off_ms;
n_frames = round(total_ms / 1000 * framerate);
onset_f  = round(onset_ms / 1000 * framerate);
offset_f = round((onset_ms + duration_ms) / 1000 * framerate);

x_deg = ((1:N) - (N+1)/2) / (N/dov);
y_deg = ((1:N) - (N+1)/2) / (N/dov);

%% BUILD STIMULUS AND FIT THE READOUT ANNULUS
I_fwd    = make_CentralDrift_stimulus_2D(N, contrast, meanL, n_segments);
I_fliplr = fliplr(I_fwd);
I_flipud = flipud(I_fwd);

[Xg, Yg] = meshgrid(1:N, 1:N);
r_map = sqrt((Xg-cx_wheel).^2 + (Yg-cy_wheel).^2);
disc_mask = abs(I_fwd - meanL) > 2;
r_outer_ro = prctile(r_map(disc_mask), 98);

fprintf('\n==== RESULTS CENTRAL DRIFT (2D) ====\n');
fprintf('N=%d, n_segments=%d, sf0_factor=%d (period=100px), sigmaW=%.2f, Nw=%dms, profile=%s\n', ...
    N, n_segments, model.sf0_factor, model.sigmaW, model.Nw, sigma_profile);
fprintf('annulus: [%.0f, %.0f) px around canvas centre\n\n', r_inner_ro, r_outer_ro);
fprintf('  %-16s %-16s\n', 'Condition', 'proj_tan');
fprintf('  %s\n', repmat('-', 1, 34));

cond_names = {'forward', 'invertido'};
n_cond = numel(cond_names);
proj_tan_result = zeros(1, n_cond);

figure(1)
set(gcf, 'Position', [50, 100, 320*n_cond, 340]);

%% RUN
for ci = 1:n_cond
    cname = cond_names{ci};

    switch cname
        case 'forward'
            TCB = single(meanL * ones(N, N, n_frames));
            for t = onset_f : min(offset_f, n_frames)
                TCB(:,:,t) = single(I_fwd);
            end
            [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
            proj_raw = make_tangential_field(O_H, O_V, cx_wheel, cy_wheel);

        case 'invertido'
            TCB_lr = single(meanL * ones(N, N, n_frames));
            TCB_ud = single(meanL * ones(N, N, n_frames));
            for t = onset_f : min(offset_f, n_frames)
                TCB_lr(:,:,t) = single(I_fliplr);
                TCB_ud(:,:,t) = single(I_flipud);
            end
            [O_H_lr, O_V_lr] = stINRF3_compute(TCB_lr, dov, model, framerate, N);
            [O_H_ud, O_V_ud] = stINRF3_compute(TCB_ud, dov, model, framerate, N);
            proj_lr = make_tangential_field(O_H_lr, O_V_lr, cx_wheel, cy_wheel);
            proj_ud = make_tangential_field(O_H_ud, O_V_ud, cx_wheel, cy_wheel);
            proj_raw = (proj_lr + proj_ud) / 2;
            TCB = TCB_lr;   % for display/video only
    end

    proj_tan_result(ci) = smooth_and_score_proj(proj_raw, ...
        cx_wheel, cy_wheel, cx_wheel, cy_wheel, r_inner_ro, r_outer_ro, sigma_profile);

    fprintf('  %-16s %+.5e\n', cname, proj_tan_result(ci));

    subplot(1, n_cond, ci)
    imagesc(x_deg, y_deg, TCB(:,:,round((onset_f+offset_f)/2)))
    colormap(gray)
    axis image
    set(gca, 'YDir', 'normal')
    xlabel('x (deg)')
    ylabel('y (deg)')
    title(cname, 'Interpreter', 'none')
    pause(0.01)

    video_filename = sprintf('CentralDrift_%s.avi', cname);
    v = VideoWriter(video_filename, 'Motion JPEG AVI');
    v.FrameRate = 60;
    v.Quality   = 100;
    open(v);
    TCB_normalized = TCB - min(TCB(:));
    denom = max(TCB_normalized(:));
    if denom == 0, denom = 1; end
    TCB_uint8 = uint8(255 * (TCB_normalized / denom));
    for tt = 1:n_frames
        writeVideo(v, TCB_uint8(:,:,tt));
    end
    close(v);
    fprintf('Video saved: %s\n', video_filename);
end

fprintf('\n');

%% PLOT
fig = figure(2);
set(fig, 'Color', 'w', 'Position', [100 100 700 500]);
bar(proj_tan_result, 'FaceColor', [0.35 0.55 0.85], 'EdgeColor', [0 0 0], 'LineWidth', 1.2);
ax = gca;
ax.XTick = 1:n_cond;
ax.XTickLabel = cond_names;
ax.TickLabelInterpreter = 'none';
ylabel('proj\_tan (a.u., + = CW)', 'FontSize', 12);
title('Central Drift — tangential projection by condition', ...
    'FontSize', 12, 'FontWeight', 'normal');
grid on; box on;
set(gca, 'FontSize', 10, 'LineWidth', 1.0, 'TickDir', 'out');

%% SAVE
outdir = fullfile(pwd, 'results_CentralDrift_2D');
if ~exist(outdir, 'dir'), mkdir(outdir); end

stamp = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(outdir, sprintf('CentralDrift_2D_%s.mat', stamp)), ...
    'proj_tan_result', 'cond_names', 'N', 'n_segments', 'dov', 'framerate', ...
    'template_path', 'sigma_profile', 'r_inner_ro', 'r_outer_ro', ...
    'cx_wheel', 'cy_wheel', 'model');

saveas(fig, fullfile(outdir, sprintf('CentralDrift_2D_%s.png', stamp)));
fprintf('Saved results to: %s\n', outdir);
