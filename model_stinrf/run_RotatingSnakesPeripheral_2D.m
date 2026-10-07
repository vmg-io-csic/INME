%% run_RotatingSnakesPeripheral_2D.m
%
% Runs the Rotating Snakes illusion (Kitaoka & Ashida, 2003) using the
% wheel GENERATED FROM A FORMULA fitted to the reference Kitaoka images
% (make_RotatingSnakesWheel_stimulus_2D.m, same generator as
% run_RotatingSnakesWheel_2D.m -- see that script's header for details),
% placed off-centre so the wheel falls in PERIPHERAL vision while the
% observer (and the model) "fixate" the canvas centre (the fovea).
% Validated design for the peripheral configuration
% Sec. 1.7). The wheel keeps the EXACT SAME pixel size (N_wheel x
% N_wheel) as the centred experiment; only the canvas grows and the
% wheel is shifted within it.
%
% Four conditions (same definitions as run_RotatingSnakesWheel_2D.m):
%   forward   — the pattern as-is, placed off-centre. Looks CW. Expected:
%               positive.
%   invertido — average of the fliplr- and flipud-reflected WHEEL
%               (reflected in place, then re-placed at the same
%               off-centre position) tangential-projection fields,
%               averaged before smoothing. Looks CCW. Expected: negative.
%   symmetric — same ring layout with a symmetric luminance profile,
%               placed off-centre, no fliplr/flipud averaging.
%               Zero-control. Expected: approx 0.
%   baseline  — uniform grey canvas throughout. Expected: approx 0.
%
% Read-out: same global tangential-projection metric as
% run_RotatingSnakesWheel_2D.m (make_tangential_field.m +
% smooth_and_score_proj.m; positive = CW, negative = CCW -- see that
% file), but with the smoothing sigma profile centred on the FOVEA
% (canvas centre) while the annulus mask stays centred on the WHEEL --
% this is where centred vs peripheral differ mechanically (report
% Sec. 1.4: the exponential profile decays from the fovea, so it is
% already small/local by the time it reaches the eccentric wheel).
%
% ECCENTRICITY: the the eccentricity (~75 px, ~0.44 deg) puts
% the wheel (radius 147 px, ~0.86 deg) close enough to the fovea that
% profile D's sigma is still ~80 px (bigger than a whole ring, ~30-45 px)
% over part of the wheel -- that heavy, sharply non-uniform blur across
% the wheel's own extent injects a spurious net tangential bias of
% geometric origin, comparable to or larger than the true signal, which
% can flip the metric's sign (verified empirically: raw unsmoothed signal
% and wheel-centred-smoothing both agree with the centred condition's
% sign; only fovea-centred smoothing at this small eccentricity flips
% it -- confirmed not a discretisation artifact by sweeping
% variable_gauss's n_levels from 12 to 400 with no material change).
% Sweeping the fovea-to-wheel separation shows the smoothed value
% converges to the raw (correct-sign) value from ~225 px onward. This
% script therefore uses ecc = 250 px (N=800) so the wheel sits fully
% outside the near-foveal high-sigma zone -- a genuinely peripheral
% placement, not an arbitrary fix to force a sign.
%
% Model parameters (report Sec. 1.1, 1.7 -- sf0_factor generalised to
% N/100 for this bigger canvas): sf0_factor = 8 for N=800 (keeps the
% Gabor period at 100 px, exactly as in the centred N=300 case where
% sf0_factor=3); sigmaMu = 0.12 deg; sigmaW = 0.22 deg (validated
% optimum); Nw = 75 ms.
%
% Expected results: forward positive (CW), invertido negative (CCW),
% roughly equal magnitude, baseline approx 0, and stronger in magnitude
% than the centred case (report: ~3.5x stronger at ecc=75px, consistent
% with the illusion being perceived more strongly in peripheral vision;
% at ecc=250px the raw signal was measured ~2.9x stronger than centred).
% n_units=23 (see rotating_snakes_params.m), which avoids the degenerate
% 90 degree symmetry of a unit count that is a multiple of four.
%
% For each condition, the stimulus is displayed and saved as a video
% (with a marker at the fovea), and a final bar plot of the tangential-
% projection scalar is saved to results_RotatingSnakesPeripheral_2D/.
%
% Inputs
%   params_illusion, params_symmetric (char) paths to the fitted ring
%             the wheel parameters (rotating_snakes_params.m)
%   N_wheel     (int)    wheel size [px], square — SAME as the centred
%                        experiment's N
%   N           (int)    canvas size [px], square (bigger than N_wheel,
%                        to leave room for off-centre placement)
%   dov         (double) field of view [dva]
%   framerate   (int)    temporal sampling [fps]
%   onset_ms, duration_ms, post_off_ms (int) temporal paradigm [ms]
%   meanL       (double) background luminance [a.u.]
%   sigma_profile (char) smoothing profile, see smooth_and_score_proj.m
%   cx_wheel, cy_wheel (double) off-centre wheel position [px]
%   cx_fovea, cy_fovea (double) fixation point (canvas centre) [px]
%

clear; clc; close all;

params_illusion  = 'illusion';
params_symmetric = 'symmetric';

%% PARAMETERS
N_wheel     = 300;      % wheel size [px] -- identical to the centred experiment
N           = 800;      % canvas size [px] -- big enough for a genuinely
                         % peripheral placement (see header note on
                         % eccentricity and the near-foveal smoothing artifact)
dov         = N * (1.5 / 128);
framerate   = 200;
onset_ms    = 150;
duration_ms = 150;
post_off_ms = 200;
meanL       = 50;

sigma_profile = 'D';    % validated optimum

% Fovea = canvas centre
cx_fovea = N / 2;
cy_fovea = N / 2;

% Wheel centre: OFF-CENTRE (peripheral), placed 250 px from the fovea so
% the wheel (radius 147 px) sits fully outside profile D's near-foveal
% high-sigma zone (see header note) -- leaves margin to stay inside the
% N=800 canvas.
cx_wheel = cx_fovea - 245;
cy_wheel = cy_fovea;
ecc_px = sqrt((cx_wheel - cx_fovea)^2 + (cy_wheel - cy_fovea)^2);

%% MODEL PARAMETERS 
cfg   = stINRF3_config('rotating_snakes');
model = cfg.model;
model.sf0    = [];
model = rmfield(model, 'sf0');
model.sf0_factor = round(N / 100);   % = 6 for N=600, keeps period at 100 px
model.sigmaMu    = 0.12;
model.sigmaW     = 0.22;             % validated optimum
model.Nw         = 75;               % adjusted from the 15 ms default

%% DERIVED
total_ms = onset_ms + duration_ms + post_off_ms;
n_frames = round(total_ms / 1000 * framerate);
onset_f  = round(onset_ms / 1000 * framerate);
offset_f = round((onset_ms + duration_ms) / 1000 * framerate);

%% BUILD WHEELS, PLACE ON THE BIGGER CANVAS, FIT THE READOUT ANNULUS
[~, I_wheel] = make_RotatingSnakesWheel_stimulus_2D(params_illusion, N_wheel, 1, 1, 1, meanL);
[~, I_wheel_sym] = make_RotatingSnakesWheel_stimulus_2D(params_symmetric, N_wheel, 1, 1, 1, meanL);
I_wheel_lr = fliplr(I_wheel);
I_wheel_ud = flipud(I_wheel);

r_start = round(cy_wheel - N_wheel/2);
c_start = round(cx_wheel - N_wheel/2);
r_end   = r_start + N_wheel - 1;
c_end   = c_start + N_wheel - 1;
if r_start < 1 || c_start < 1 || r_end > N || c_end > N
    error('Wheel does not fit in canvas. Increase N or move cx_wheel/cy_wheel.');
end

place_on_canvas = @(wheel_img) local_place(wheel_img, N, meanL, r_start, r_end, c_start, c_end);
I_fwd    = place_on_canvas(I_wheel);
I_sym    = place_on_canvas(I_wheel_sym);
I_fliplr = place_on_canvas(I_wheel_lr);
I_flipud = place_on_canvas(I_wheel_ud);

[Xg, Yg] = meshgrid(1:N, 1:N);
r_map = sqrt((Xg-cx_wheel).^2 + (Yg-cy_wheel).^2);
disc_mask = abs(I_fwd - meanL) > 2;
r_outer_ro = prctile(r_map(disc_mask), 98);
r_inner_ro = 0;   % include the full disc -- excluding the centre changes
                   % the metric by ~1% (checked empirically), not worth
                   % the extra parameter

% Pre-calculate physical units for axes labeling
x_deg = ((1:N) - (N+1)/2) / (N/dov);
y_deg = ((1:N) - (N+1)/2) / (N/dov);

fprintf('\n==== RESULTS ROTATING SNAKES PERIPHERAL (2D, parametric stimulus) ====\n');
fprintf('N=%d, sf0_factor=%d (period=100px), sigmaW=%.2f, Nw=%dms, profile=%s\n', ...
    N, model.sf0_factor, model.sigmaW, model.Nw, sigma_profile);
fprintf('wheel=(%d,%d), fovea=(%d,%d), ecc=%.0f px, annulus=[%.0f,%.0f) px\n\n', ...
    cx_wheel, cy_wheel, cx_fovea, cy_fovea, ecc_px, r_inner_ro, r_outer_ro);
fprintf('  %-16s %-16s\n', 'Condition', 'proj_tan');
fprintf('  %s\n', repmat('-', 1, 34));

cond_names = {'forward', 'invertido', 'symmetric', 'baseline'};
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
            proj_raw = (proj_lr + proj_ud) / 2;   % averaged BEFORE smoothing
            TCB = TCB_lr;   % for display/video only

        case 'symmetric'
            TCB = single(meanL * ones(N, N, n_frames));
            for t = onset_f : min(offset_f, n_frames)
                TCB(:,:,t) = single(I_sym);
            end
            [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
            proj_raw = make_tangential_field(O_H, O_V, cx_wheel, cy_wheel);

        case 'baseline'
            TCB = single(meanL * ones(N, N, n_frames));
            [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
            proj_raw = make_tangential_field(O_H, O_V, cx_wheel, cy_wheel);
    end

    proj_tan_result(ci) = smooth_and_score_proj(proj_raw, ...
        cx_wheel, cy_wheel, cx_fovea, cy_fovea, r_inner_ro, r_outer_ro, sigma_profile);

    fprintf('  %-16s %+.5e\n', cname, proj_tan_result(ci));

    subplot(1, n_cond, ci)
    imagesc(x_deg, y_deg, TCB(:,:,round((onset_f+offset_f)/2)))
    colormap(gray)
    axis image
    set(gca, 'YDir', 'normal')
    hold on
    plot((cx_fovea-(N+1)/2)/(N/dov), (cy_fovea-(N+1)/2)/(N/dov), 'g+', 'MarkerSize', 10, 'LineWidth', 1.5)
    hold off
    xlabel('x (deg)')
    ylabel('y (deg)')
    title(cname, 'Interpreter', 'none')
    pause(0.01)

    video_filename = sprintf('RotSnakesPeri_%s.avi', cname);
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
title('Rotating Snakes Peripheral — tangential projection by condition', ...
    'FontSize', 12, 'FontWeight', 'normal');
grid on; box on;
set(gca, 'FontSize', 10, 'LineWidth', 1.0, 'TickDir', 'out');

%% SAVE
outdir = fullfile(pwd, 'results_RotatingSnakesPeripheral_2D');
if ~exist(outdir, 'dir'), mkdir(outdir); end

stamp = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(outdir, sprintf('RotSnakesPeri_2D_%s.mat', stamp)), ...
    'proj_tan_result', 'cond_names', 'N', 'dov', 'framerate', 'sigma_profile', ...
    'params_illusion', 'params_symmetric', ...
    'cx_wheel', 'cy_wheel', 'cx_fovea', 'cy_fovea', 'ecc_px', ...
    'r_inner_ro', 'r_outer_ro', 'model');

saveas(fig, fullfile(outdir, sprintf('RotSnakesPeri_2D_%s.png', stamp)));
fprintf('Saved results to: %s\n', outdir);

%% LOCAL FUNCTIONS
function canvas = local_place(wheel_img, N, meanL, r_start, r_end, c_start, c_end)
    canvas = meanL * ones(N, N);
    canvas(r_start:r_end, c_start:c_end) = wheel_img;
end
