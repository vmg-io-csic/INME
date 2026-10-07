%% run_RotatingSnakesWheel_2D.m
%
% Four conditions:
%   forward    the wheel looks CW.
%   symmetric  the symmetric-control wheel
%   baseline   uniform grey canvas throughout. Expected: approx 0.
%

% Inputs
%   params_illusion, params_symmetric (char) paths to the fitted ring
%             the wheel parameters (rotating_snakes_params.m)
%   N              (int)    canvas size [px], square
%   dov            (double) field of view [dva]
%   framerate      (int)    temporal sampling [fps]
%   onset_ms       (int)    pattern onset time [ms]
%   duration_ms    (int)    pattern visible duration [ms]
%   post_off_ms    (int)    frames after offset [ms]
%   meanL          (double) background luminance [a.u.]
%   sigma_profile  (char)   smoothing profile
%   r_inner_ro, r_outer_ro (double) read-out  radii [px]
%

clear; clc; close all;

params_illusion  = 'illusion';
params_symmetric = 'symmetric';

%% STIMULUS PARAMETERS
N           = 300;      % canvas size [px]
dov         = N * (1.5 / 128);
framerate   = 200;
onset_ms    = 150;
duration_ms = 150;
post_off_ms = 200;
meanL       = 50;

sigma_profile = 'D';    % validated optimum (validated optimum)

%% MODEL PARAMETERS (see the model parameter table)
cfg   = stINRF3_config('rotating_snakes');
model = cfg.model;
model.sf0    = [];
model = rmfield(model, 'sf0');
model.sf0_factor = round(N / 100);   % keeps the Gabor period at 100 px
model.sigmaMu    = 0.12;
model.sigmaW     = 0.22;             % validated optimum
model.Nw         = 75;               % adjusted from the 15 ms default

%% READOUT
r_inner_ro = 0;   % include the full disc -- excluding the centre changes
                   % the metric by ~1% (checked empirically), not worth
                   % the extra parameter
cx_wheel   = (N+1)/2;
cy_wheel   = (N+1)/2;

%% DERIVED
total_ms = onset_ms + duration_ms + post_off_ms;
n_frames = round(total_ms / 1000 * framerate);
onset_f  = round(onset_ms / 1000 * framerate);
offset_f = round((onset_ms + duration_ms) / 1000 * framerate);

% Pre-calculate physical units for axes labeling
x_deg = ((1:N) - (N+1)/2) / (N/dov);
y_deg = ((1:N) - (N+1)/2) / (N/dov);

[~, I_fwd] = make_RotatingSnakesWheel_stimulus_2D(params_illusion, N, 1, 1, 1, meanL);
[~, I_sym] = make_RotatingSnakesWheel_stimulus_2D(params_symmetric, N, 1, 1, 1, meanL);
I_fliplr = fliplr(I_fwd);
I_flipud = flipud(I_fwd);

[Xg, Yg] = meshgrid(1:N, 1:N);
r_map = sqrt((Xg-cx_wheel).^2 + (Yg-cy_wheel).^2);
disc_mask = abs(I_fwd - meanL) > 2;
r_outer_ro = prctile(r_map(disc_mask), 98);

fprintf('\n==== RESULTS ROTATING SNAKES WHEEL (2D) ====\n');
fprintf('N=%d, sf0_factor=%d (period=100px), sigmaW=%.2f, Nw=%dms, profile=%s\n', ...
    N, model.sf0_factor, model.sigmaW, model.Nw, sigma_profile);
fprintf('annulus: [%.0f, %.0f) px around wheel centre\n\n', r_inner_ro, r_outer_ro);
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

    video_filename = sprintf('RotSnakesWheel_%s.avi', cname);
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
title('Rotating Snakes Wheel — tangential projection by condition', ...
    'FontSize', 12, 'FontWeight', 'normal');
grid on; box on;
set(gca, 'FontSize', 10, 'LineWidth', 1.0, 'TickDir', 'out');

%% SAVE
outdir = fullfile(pwd, 'results_RotatingSnakesWheel_2D');
if ~exist(outdir, 'dir'), mkdir(outdir); end

stamp = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(outdir, sprintf('RotSnakesWheel_2D_%s.mat', stamp)), ...
    'proj_tan_result', 'cond_names', 'N', 'dov', 'framerate', ...
    'params_illusion', 'params_symmetric', 'sigma_profile', 'r_inner_ro', 'r_outer_ro', ...
    'cx_wheel', 'cy_wheel', 'model');

saveas(fig, fullfile(outdir, sprintf('RotSnakesWheel_2D_%s.png', stamp)));
fprintf('Saved results to: %s\n', outdir);
