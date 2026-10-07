%% analyse_interaction_threshold.m
%
% Locates the envelope width at which the CSr condition crosses chance,
% reading the tensor saved by run_Interaction_1D. Nothing is recomputed,
% so this runs in a second on the login node.
%
% Two independent estimates are reported:
%   - the width at which the percentage of correctly signed trials crosses
%     50%, interpolated linearly in log(sigma_t)
%   - the width at which the mean response changes sign
% They should agree; a large gap means the distribution across phases is
% skewed there and the percentage is the less reliable estimator.
%
% Inputs
%   mat_path (char) results file to read; empty picks the most recent
%   cond_of_interest (char) condition whose threshold is plotted in detail
%

clear; clc; close all;

mat_path = '';                    % empty = most recent results/interaction_*.mat
cond_of_interest = 'CSr';

%% LOAD
if isempty(mat_path)
    here = fileparts(mfilename('fullpath'));
    if isempty(here), here = pwd; end
    d = dir(fullfile(here, 'results', 'interaction_*.mat'));
    if isempty(d)
        error('No results/interaction_*.mat found. Run run_Interaction_1D first.');
    end
    [~, newest] = max([d.datenum]);
    mat_path = fullfile(d(newest).folder, d(newest).name);
end
fprintf('reading %s\n', mat_path);
S = load(mat_path);

val        = S.val;                       % n_cond x n_sigma x n_trials
sigma_list = S.sigma_t_list(:).';
names      = S.cond_names;
direction  = S.direction;
[n_cond, n_sigma, n_trials] = size(val);

fprintf('%d conditions x %d envelope widths x %d trials\n', n_cond, n_sigma, n_trials);
fprintf('sigma_t from %g to %g ms\n\n', min(sigma_list), max(sigma_list));

%% METRICS
pct  = 100 * mean(val * direction > 0, 3);
mu   = mean(val, 3);
sem  = std(val, 0, 3) / sqrt(n_trials);

%% TABLES
print_table('Sign matches stimulus direction (%)', names, sigma_list, pct,  '%11.0f%%');
print_table('Mean response (a.u.)',                names, sigma_list, mu,   '%12.4f');
print_table('SEM of the response (a.u.)',          names, sigma_list, sem,  '%12.4f');

%% THRESHOLDS
fprintf('\n--- Chance-level crossing (50%% correct) ---\n');
thr_pct = nan(n_cond,1);
for ic = 1:n_cond
    thr_pct(ic) = crossing(sigma_list, pct(ic,:), 50);
    report(names{ic}, thr_pct(ic), sigma_list, pct(ic,:), 50);
end

fprintf('\n--- Sign reversal of the mean response ---\n');
thr_mu = nan(n_cond,1);
for ic = 1:n_cond
    thr_mu(ic) = crossing(sigma_list, mu(ic,:)*direction, 0);
    report(names{ic}, thr_mu(ic), sigma_list, mu(ic,:)*direction, 0);
end

%% FIGURE: all conditions, then the one of interest on its own scale
ic = find(strcmp(names, cond_of_interest), 1);
if isempty(ic), ic = 1; end

figure('Color','w','Position',[80 80 1000 420]);

subplot(1,2,1); hold on;
markers = {'o','s','d','^','v','>','<'};
for c = 1:n_cond
    errorbar(sigma_list, mu(c,:), sem(c,:), ['-' markers{mod(c-1,numel(markers))+1}], ...
             'LineWidth', 1.3, 'MarkerSize', 6, 'MarkerFaceColor', 'w');
end
yline(0,'k--','LineWidth',1);
xlabel('temporal envelope \sigma_t (ms)'); ylabel('mean response (a.u.)');
title('all conditions'); legend(names, 'Location','best'); grid on; box on;

subplot(1,2,2); hold on;
fill([sigma_list fliplr(sigma_list)], ...
     [mu(ic,:)+sem(ic,:) fliplr(mu(ic,:)-sem(ic,:))], [0.2 0.4 0.75], ...
     'FaceAlpha', 0.25, 'EdgeColor','none');
plot(sigma_list, mu(ic,:), '-o', 'Color', [0.2 0.4 0.75], 'LineWidth', 1.6, ...
     'MarkerFaceColor','w');
yline(0,'k--','LineWidth',1);
if ~isnan(thr_mu(ic))
    xline(thr_mu(ic), 'r-', 'LineWidth', 1.4);
    text(thr_mu(ic), max(mu(ic,:)), sprintf('  %.1f ms', thr_mu(ic)), ...
         'Color','r', 'VerticalAlignment','top');
end
xlabel('temporal envelope \sigma_t (ms)'); ylabel('mean response (a.u.)');
title(sprintf('%s, own scale', names{ic})); grid on; box on;

fprintf('\n');

%% LOCAL FUNCTIONS
function print_table(header, names, sigma_list, M, fmt)
% PRINT_TABLE  One metric as conditions x envelope widths.
    fprintf('\n--- %s ---\n', header);
    fprintf('  %-6s', 'Cond');
    for i = 1:numel(sigma_list), fprintf('%12s', sprintf('%g ms', sigma_list(i))); end
    fprintf('\n  %s\n', repmat('-', 1, 6 + 12*numel(sigma_list)));
    for c = 1:numel(names)
        fprintf('  %-6s', names{c});
        fprintf(fmt, M(c,:));
        fprintf('\n');
    end
end

function x = crossing(xs, ys, level)
% CROSSING  Last upward crossing of level, interpolated in log(x).
    k = find(ys(1:end-1) < level & ys(2:end) >= level, 1, 'last');
    if isempty(k)
        x = NaN;
        return;
    end
    l1 = log(xs(k));  l2 = log(xs(k+1));
    x  = exp(l1 + (level - ys(k)) * (l2 - l1) / (ys(k+1) - ys(k)));
end

function report(name, thr, xs, ys, level)
% REPORT  One line describing where a curve crosses a level.
    if ~isnan(thr)
        k = find(ys(1:end-1) < level & ys(2:end) >= level, 1, 'last');
        fprintf('  %-6s %8.2f ms   (between %g and %g)\n', name, thr, xs(k), xs(k+1));
    elseif all(ys >= level)
        fprintf('  %-6s always at or above the level across the sweep\n', name);
    elseif all(ys < level)
        fprintf('  %-6s always below the level across the sweep\n', name);
    else
        fprintf('  %-6s no clean upward crossing (non-monotonic)\n', name);
    end
end
