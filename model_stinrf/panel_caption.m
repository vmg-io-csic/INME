function panel_caption(ax, name, d, c)
% PANEL_CAPTION  Two-line caption under an image panel: the condition name
%   in dark grey and, below it, the direction the stimulus moves, in the
%   colour of the frame.
%
%   Drawn with annotation in figure coordinates rather than with text in
%   axes coordinates, because the automatic x label takes an unpredictable
%   amount of room and a caption placed relative to the axes ends up either
%   on top of it or off the canvas.
%
% Inputs
%   ax (axes) panel the caption belongs to
%   name (char) condition name
%   d (double) direction the stimulus moves, positive rightward
%   c (struct) palette from stinrf_fig_style
%
% Outputs
%   nothing; draws into the parent figure

    fig = ancestor(ax, 'figure');
    p   = get(ax, 'Position');

    annotation(fig, 'textbox', [p(1), p(2)-0.155, p(3), 0.055], ...
        'String', name, 'FontName','Arial', 'FontSize', 10, ...
        'Color', c.text, 'EdgeColor','none', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');

    if d >= 0
        word = 'Rightward';  col = c.right;
    else
        word = 'Leftward';   col = c.left;
    end
    annotation(fig, 'textbox', [p(1), p(2)-0.215, p(3), 0.055], ...
        'String', word, 'FontName','Arial', 'FontSize', 9, ...
        'FontWeight','bold', 'Color', col, 'EdgeColor','none', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
