function plot_runs_bar(ax, pos, runs, col, opts)
% PLOT_RUNS_BAR  Mean as a bar, every run as a dot on top of it.
%
%   The house way of showing a multi-run condition: a semi-transparent bar
%   for the mean, one translucent dot per run jittered horizontally so they
%   do not pile up, and a black error bar for the standard deviation. The
%   dots are translucent so that they read over both the bar and the error
%   bar without hiding either.
%
% Inputs
%   ax (axes) target axes
%   pos (double) x position of the bar
%   runs (1xR double) one value per run
%   col (1x3 double) bar colour
%   opts (struct, optional) fields:
%       width (double)  bar width, default 0.55
%       jitter (double) half-width of the dot scatter, default 0.13
%       dot (double)    dot area, default 16
%       alpha (double)  dot opacity, default 0.45
%       err (char)      'sd' (default) or 'sem'
%
% Outputs
%   nothing; draws into ax

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'width'),  opts.width  = 0.55; end
    if ~isfield(opts,'jitter'), opts.jitter = 0.13; end
    if ~isfield(opts,'dot'),    opts.dot    = 16;   end
    if ~isfield(opts,'alpha'),  opts.alpha  = 0.45; end
    if ~isfield(opts,'err'),    opts.err    = 'sd'; end

    runs = runs(:).';
    m = mean(runs);
    if strcmp(opts.err, 'sem')
        e = std(runs) / sqrt(numel(runs));
    else
        e = std(runs);
    end

    hold(ax, 'on');
    bar(ax, pos, m, opts.width, 'FaceColor', col, 'FaceAlpha', 0.55, ...
        'EdgeColor', col, 'LineWidth', 1.2);

    jit = (rand(1, numel(runs)) - 0.5) * 2 * opts.jitter;
    scatter(ax, pos + jit, runs, opts.dot, 'filled', ...
            'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerFaceAlpha', opts.alpha, 'MarkerEdgeColor', 'none');

    errorbar(ax, pos, m, e, 'k', 'LineStyle', 'none', 'LineWidth', 1.3, ...
             'CapSize', 7);
end
