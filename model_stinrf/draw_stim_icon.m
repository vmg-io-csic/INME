function draw_stim_icon(ax, comp_dirs, pattern_dir, cols, cycles, amp)
% DRAW_STIM_ICON  Circular stimulus icon in the style of Movshon plaid figures.
%
%   A light grey disc carrying square-wave components, with one arrow per
%   component and, when given, a longer arrow for the pattern direction.
%   The wave vector of each component is aligned with its drift direction,
%   so the bars come out perpendicular to their own arrow.
%
% Inputs
%   ax (axes) target axes
%   comp_dirs (1xK double) component drift directions [deg], K = 1 or 2
%   pattern_dir (double) pattern direction [deg], [] to omit
%   cols (struct, optional) fields comp1, comp2, pattern
%   cycles (double, optional) periods across the disc; a scalar applies to
%       every component, a vector gives one value per component so that two
%       components of different spatial frequency can be shown superimposed.
%       Values above 8 are clipped, since finer bars alias at icon size
%   amp (double, optional) modulation depth about mid grey, 0 to 1, so that
%       a series of icons can show how contrast varies. Default 1
%
% Outputs
%   nothing; draws into ax

    if nargin < 4 || isempty(cols)
        cols.comp1   = [0.00 0.62 0.45];
        cols.comp2   = [0.90 0.62 0.00];
        cols.pattern = [0.00 0.45 0.70];
    end
    if nargin < 3, pattern_dir = []; end

    if nargin < 6 || isempty(amp), amp = 1; end
    if nargin < 5 || isempty(cycles), cycles = 2.2; end
    if isscalar(cycles), cycles = repmat(cycles, 1, numel(comp_dirs)); end
    cycles = min(cycles, 8);

    duty = 0.42;
    bg   = 0.86;          % light grey outside the bars, as in the reference

    n = 420;
    v = linspace(-1, 1, n);
    [X, Y] = meshgrid(v, v);
    inside = (X.^2 + Y.^2) <= 1;

    L = zeros(n);
    for k = 1:numel(comp_dirs)
        th = deg2rad(comp_dirs(k));
        u  = X*cos(th) + Y*sin(th);
        L  = L + double(mod(u * cycles(k), 1) >= duty);
    end
    img = 0.5*bg + amp * bg * (L/numel(comp_dirs) - 0.5);
    img(~inside) = NaN;

    cla(ax); hold(ax, 'on');
    imagesc(ax, v, v, img, 'AlphaData', ~isnan(img));
    colormap(ax, gray(256));
    set(ax, 'CLim', [0 1], 'YDir', 'normal');

    a = linspace(0, 2*pi, 400);
    plot(ax, cos(a), sin(a), 'k-', 'LineWidth', 1.0);

    cc = {cols.comp1, cols.comp2};
    for k = 1:numel(comp_dirs)
        th = deg2rad(comp_dirs(k));
        arrow(ax, 1.08*[cos(th) sin(th)], 1.62*[cos(th) sin(th)], ...
              cc{min(k,2)}, 1.8, 0.34);
    end
    if ~isempty(pattern_dir)
        th = deg2rad(pattern_dir);
        arrow(ax, 1.08*[cos(th) sin(th)], 2.05*[cos(th) sin(th)], ...
              cols.pattern, 2.6, 0.42);
    end

    axis(ax, 'equal');
    axis(ax, [-2.35 2.35 -2.35 2.35]);
    axis(ax, 'off');
    hold(ax, 'off');
end

function arrow(ax, p0, p1, col, lw, head)
% ARROW  Line with a filled triangular head, in data coordinates.
%   head is the length of the head; it is made generous on purpose, because
%   at icon scale a small head disappears.
    d = (p1 - p0) / max(norm(p1 - p0), eps);
    stem_end = p1 - 0.85*head*d;
    plot(ax, [p0(1) stem_end(1)], [p0(2) stem_end(2)], '-', ...
         'Color', col, 'LineWidth', lw);
    n = [-d(2) d(1)];
    base = p1 - head*d;
    patch(ax, 'XData', [p1(1), base(1)+0.5*head*n(1), base(1)-0.5*head*n(1)], ...
              'YData', [p1(2), base(2)+0.5*head*n(2), base(2)-0.5*head*n(2)], ...
              'FaceColor', col, 'EdgeColor', col, 'LineWidth', 0.5);
end
