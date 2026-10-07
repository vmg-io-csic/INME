%% run_Interaction_examples_1D.m
%
% Figure for the cross-scale interaction: the CSr condition, in which the
% 3 cpd component drifts and the 1 cpd component is static, reverses the
% reported direction when the temporal envelope is short.
%
% Three envelope widths are shown: short, where the direction is inverted;
% the width at which performance is at chance; and long, where the
% direction is correct. Each has its space-time panel, and the bar chart
% carries one dot per run plus the standard deviation.
%
% Sign convention: the stimulus drifts leftward, so a negative response is
% the correct direction and a positive one is a reversal.
%
% Inputs
%   cond_idx (int) condition to plot, 4 = CSr
%   sigma_examples (1x3 double) envelope widths [ms], short, chance, long
%   labels (1x3 cell) names shown above each panel
%   n_trials (int) random-phase repetitions
%   space_win_deg (double) half-width of the space axis in the panels [deg]
%
% Outputs (in results/):
%   Fig_Interaction.png (600 dpi) and .pdf
%   interaction_examples_<stamp>.mat
%

clear; clc; close all;

cfg   = stINRF3_config('DH');
model = cfg.model;

cond_idx       = 4;
sigma_examples = [10 89.2 200];
labels         = {'', '', ''};        % the width itself is the label
space_win_deg  = 2;

N            = 400;
dov          = 4.0;
framerate    = 1000;
duration_ms  = 500;
meanL        = 100;
sfm_Low_cpd  = 1;
sfm_High_cpd = 3;
direction    = -1;
tf_hz        = 12;
v_Low  = tf_hz / sfm_Low_cpd;
v_High = tf_hz / sfm_High_cpd;

n_trials = 25;
rng_seed = 42;
sigma_c  = 57.74;

out_dir = stinrf_results_dir(mfilename);
if isempty(out_dir), out_dir = pwd; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

n_frames  = round(duration_ms/1000*framerate);
t_ms      = (0:n_frames-1)/framerate*1000;
x_deg     = (0:N-1) * dov / N - dov/2;
j0        = round(N/2);
j_sensors = (j0-99):(j0+100);
nE        = numel(sigma_examples);
rng(rng_seed);
ph_high = rand(nE, n_trials) * 2*pi;
ph_low  = rand(nE, n_trials) * 2*pi;

fprintf('\n==== CROSS-SCALE INTERACTION: CSr ====\n');
fprintf('%d envelope widths x %d trials | %d fps | leftward\n', ...
        nE, n_trials, framerate);

n_jobs = nE * n_trials;
[ie_g, it_g] = ndgrid(1:nE, 1:n_trials);
ie_vec = ie_g(:);  it_vec = it_g(:);
val_flat = zeros(1, n_jobs);

parfor jj = 1:n_jobs
    e = ie_vec(jj);  it = it_vec(jj);
    mI = make_Interaction_stimulus_1D(N, n_frames, framerate, dov, ...
        cond_idx, sigma_examples(e), ph_high(e,it), ph_low(e,it), ...
        v_High, v_Low, sfm_High_cpd, sfm_Low_cpd, meanL, direction);
    O = stINRF3_compute_1D(mI, dov, model, framerate, N, j_sensors);
    val_flat(jj) = stINRF3_readout_gaussian_1D(O, j_sensors, sigma_c);
end

val = reshape(val_flat, nE, n_trials);
mu  = mean(val, 2);
sd  = std(val, 0, 2);
pct = 100 * mean(val * direction > 0, 2);

fprintf('\n  %-10s %12s %12s %10s\n', 'sigma_t', 'mean', 'SD', 'correct');
for e = 1:nE
    fprintf('  %-10s %12.5f %12.5f %9.0f%%\n', ...
            sprintf('%g ms', sigma_examples(e)), mu(e), sd(e), pct(e));
end

stim = cell(nE,1);
for e = 1:nE
    stim{e} = make_Interaction_stimulus_1D(N, n_frames, framerate, dov, ...
        cond_idx, sigma_examples(e), 0, 0, ...
        v_High, v_Low, sfm_High_cpd, sfm_Low_cpd, meanL, direction);
end

save(fullfile(out_dir, ['interaction_examples_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), ...
     'cond_idx','sigma_examples','labels','n_trials','val','mu','sd','pct', ...
     'sigma_c','model','framerate','duration_ms','stim');

%% FIGURE
c = stinrf_fig_style();
keep = abs(x_deg) <= space_win_deg;
clim = [min(cellfun(@(S) min(S(:)), stim)), max(cellfun(@(S) max(S(:)), stim))];

fig = figure('Color','w','Units','inches','Position',[1 1 10.6 4.8], ...
             'PaperPositionMode','auto');

for e = 1:nE
    ax = axes('Units','normalized','Position',[0.045 + (e-1)*0.170, 0.38, 0.135, 0.50]);
    imagesc(x_deg(keep), t_ms, stim{e}(:,keep));
    colormap(gray); set(gca,'CLim',clim,'YDir','reverse');
    xlabel('Space (deg)','FontName','Arial','FontSize',9);
    if e == 1
        ylabel('Time (ms)','FontName','Arial','FontSize',9);
    else
        set(gca,'YTickLabel',[]);
    end
    stinrf_fig_style(ax);
    xl = xlim(ax);  yl = ylim(ax);
    rectangle(ax, 'Position', [xl(1) yl(1) diff(xl) diff(yl)], ...
              'EdgeColor', c.left, 'LineWidth', 2.6);   % stimulus drifts leftward
    panel_caption(ax, sprintf('\\sigma_t = %g ms', sigma_examples(e)), ...
                  direction, c);
end

ax = axes('Units','normalized','Position',[0.640 0.38 0.325 0.50]); hold on;
for e = 1:nE
    col = dir_col_local(sign(mu(e)), c);   % green rightward, blue leftward
    plot_runs_bar(ax, e, val(e,:), col);
end
yline(0, 'k-', 'LineWidth', 0.9);
set(gca, 'XTick', 1:nE, 'XTickLabel', ...
    arrayfun(@(v) sprintf('\\sigma_t = %g ms', v), sigma_examples, ...
             'UniformOutput', false), 'TickLabelInterpreter','tex');
xlim([0.4 nE+0.6]);
ylabel('Mean response (a.u.)','FontName','Arial','FontSize',11);
set(gca,'FontSize',11);
stinrf_fig_style(ax, 'Model responses');

% exact value on every bar; the three conditions differ by orders of
% magnitude, so the numbers carry information the axis cannot show
yr = ylim(gca);  pad = 0.06*diff(yr);
for e = 1:nE
    if mu(e) >= 0, va = 'bottom'; off = pad; else, va = 'top'; off = -pad; end
    text(gca, e, mu(e)+off, sprintf('%.4f', mu(e)), 'FontName','Arial', ...
         'FontSize', 10, 'HorizontalAlignment','center', 'VerticalAlignment', va);
end
ylim(gca, yr + [-1 1]*0.16*diff(yr));

exportgraphics(fig, fullfile(out_dir,'Fig_Interaction.png'),'Resolution',600);
exportgraphics(fig, fullfile(out_dir,'Fig_Interaction.pdf'),'ContentType','vector');
fprintf('\nSaved to %s\n\n', out_dir);

function col = dir_col_local(d, c)
% DIR_COL_LOCAL  Green for rightward, blue for leftward.
    if d >= 0, col = c.right; else, col = c.left; end
end

