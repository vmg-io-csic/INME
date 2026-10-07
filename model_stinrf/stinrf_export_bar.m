function stinrf_export_bar(name, rows)
% STINRF_EXPORT_BAR  Write the raw per-phase responses beside a figure.
%
%   The mean and the SEM the figure draws are in meta.json, computed by the
%   experiment like every other figure's. This file is the measurement behind
%   them, kept so the summary can be checked, and nothing reads it to draw.
%
% Inputs
%   name (char) figure name, an existing export subfolder
%   rows (Nx3 double) phase index, forward response, reversed response

    d = getenv('STINRF_EXPORT');
    if isempty(d)
        d = fullfile(inme_output(), 'export');
    end
    d = fullfile(d, name);
    if ~exist(d, 'dir'), mkdir(d); end
    writematrix(rows, fullfile(d, 'trials.csv'));
end
