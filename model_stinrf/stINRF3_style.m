function varargout = stINRF3_style(ax, panel_letter)
% STINRF3_STYLE  House style and shared palette for the paper figures.
%
% Inputs
%   ax (axes, optional) axes to style
%   panel_letter (char, optional) bold letter drawn above the axes
%
% Outputs
%   c (struct) palette: model, inverted, ref, comp1, comp2

    if nargin == 0
        varargout{1} = palette();
        return;
    end
    set(ax, 'Box','off', 'TickDir','out', 'LineWidth',0.8, ...
            'FontName','Helvetica', 'FontSize',9, ...
            'XGrid','off', 'YGrid','off', 'Layer','top');
    if nargin > 1 && ~isempty(panel_letter)
        p = get(ax, 'Position');
        annotation(ancestor(ax,'figure'), 'textbox', ...
            [max(p(1)-0.055,0), min(p(2)+p(4)-0.01,0.94), 0.05, 0.06], ...
            'String', panel_letter, 'FontWeight','bold', 'FontSize',12, ...
            'FontName','Helvetica', 'EdgeColor','none', ...
            'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    end
    if nargout > 0, varargout{1} = palette(); end
end

function c = palette()
    c.model    = [0.15 0.20 0.55];
    c.inverted = [0.75 0.10 0.10];
    c.ref      = [0.45 0.45 0.45];
    c.comp1    = [0.10 0.55 0.25];
    c.comp2    = [0.15 0.20 0.55];
end
