function cellpress_axes(ax, xt, yt, opts)
% CELLPRESS_AXES  Offset, trimmed axes in the style used by Neuron and eLife.
%
%   Hides the native rulers and draws each axis as a single line spanning
%   only the range of the ticks, offset outwards from the data. No frame,
%   no corner, and no axis line where there are no data.
%
%   Call it last: it reads the current limits, so changing xlim or ylim
%   afterwards leaves the drawn axes behind.
%
% Inputs
%   ax (axes) target axes
%   xt (1xN double) x tick positions; [] keeps the current ones
%   yt (1xM double) y tick positions; [] keeps the current ones
%   opts (struct, optional) field pad, offset as a fraction of the range
%
% Outputs
%   nothing; modifies ax in place

    if nargin < 4 || ~isstruct(opts), opts = struct(); end
    if ~isfield(opts,'pad'), opts.pad = 0.06; end

    if ~isempty(xt), set(ax,'XTick',xt); else, xt = get(ax,'XTick'); end
    if ~isempty(yt), set(ax,'YTick',yt); else, yt = get(ax,'YTick'); end

    xr = get(ax,'XLim');  yr = get(ax,'YLim');
    dx = (xr(2)-xr(1)) * opts.pad;
    dy = (yr(2)-yr(1)) * opts.pad;
    xlo = xr(1) - dx;  ylo = yr(1) - dy;

    set(ax, 'Box','off', 'XColor','none', 'YColor','none', 'Color','none');
    hold(ax,'on');

    plot(ax, [xt(1) xt(end)], [ylo ylo], 'k-', 'LineWidth', 0.9, 'Clipping','off');
    for i = 1:numel(xt)
        plot(ax, [xt(i) xt(i)], [ylo ylo-dy*0.3], 'k-', ...
             'LineWidth', 0.9, 'Clipping','off');
        text(ax, xt(i), ylo-dy*0.5, tick_str(xt(i)), 'Clipping','off', ...
             'HorizontalAlignment','center', 'VerticalAlignment','top', ...
             'FontName','Arial', 'FontSize', 8);
    end

    plot(ax, [xlo xlo], [yt(1) yt(end)], 'k-', 'LineWidth', 0.9, 'Clipping','off');
    for i = 1:numel(yt)
        plot(ax, [xlo xlo-dx*0.3], [yt(i) yt(i)], 'k-', ...
             'LineWidth', 0.9, 'Clipping','off');
        text(ax, xlo-dx*0.45, yt(i), tick_str(yt(i)), 'Clipping','off', ...
             'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
             'FontName','Arial', 'FontSize', 8);
    end

    set(ax, 'XLim', xr, 'YLim', yr);
end

function s = tick_str(v)
% TICK_STR  Compact tick label.
    if abs(v - round(v)) < 1e-9
        s = sprintf('%d', round(v));
    else
        s = strtrim(sprintf('%.3g', v));
    end
end
