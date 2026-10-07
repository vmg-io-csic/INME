function stinrf_export(name, kind, fields, varargin)
% STINRF_EXPORT  Write one figure's data for plot_figures.py.
%
%   Every experiment ends by calling this, and every figure is drawn from
%   what it writes. One call shape for all of them is the point: when each
%   experiment exported its own way, three different conventions for the
%   error bars went unnoticed for weeks.
%
%   The export is a folder per figure holding meta.json, plus the CSVs that
%   figure needs. meta.json carries the kind, which decides the layout, and
%   whatever that kind reads.
%
% Inputs
%   name (char) figure name, becomes the export subfolder
%   kind (string) layout: "panels", "st_row", "frames_pair", "directions",
%       "curves", "masking" or "contrast"
%   fields (struct) everything else the kind needs: labels, means, SDs
%       (the whiskers the figures draw), SEMs, geometry, and so on. Written
%       into meta.json as given
%
% Name-value parameters:
%   'Matrices' (cell) space-time matrices, written as stim_1.csv onwards
%   'Curve' (Nx2 double) a measured curve, written as curve.csv
%   'OutDir' (char) export root. Defaults to the STINRF_EXPORT environment
%       variable if set, otherwise an export folder beside this file, so a
%       clone of the repository writes inside itself and needs no setting up
%
% Outputs
%   nothing; writes the folder

    p = inputParser;
    p.addParameter('Matrices', {});
    p.addParameter('Curve', []);
    p.addParameter('OutDir', default_export_dir());
    p.parse(varargin{:});
    o = p.Results;

    % A spread at the level of floating-point noise means the experiment
    % randomised nothing: every run gave the same answer. It is written as an
    % exact zero, which the figure draws as a zero-length whisker on the bar.
    % Catching it here rather than in each experiment is what keeps the
    % convention the same everywhere. The SD and the SEM are spreads of the
    % mean response; correct_sd and correct_sem of the share correct.
    pairs = {'sem', 'mean'; 'sd', 'mean'; ...
             'correct_sd', 'correct'; 'correct_sem', 'correct'};
    for k = 1:size(pairs, 1)
        sf = pairs{k, 1};  mf = pairs{k, 2};
        if ~(isfield(fields, sf) && isfield(fields, mf)), continue; end
        s = fields.(sf)(:).';
        s(~isfinite(s)) = 0;
        if max(s) <= 1e-9 * max(abs(fields.(mf)(:).')), s(:) = 0; end
        fields.(sf) = s;
        assert(numel(fields.(sf)) == numel(fields.(mf)), ...
               'stinrf_export: %d values in %s but %d in %s', ...
               numel(fields.(mf)), mf, numel(fields.(sf)), sf);
    end
    if isfield(fields, 'labels') && isfield(fields, 'mean')
        assert(numel(fields.labels) == numel(fields.mean), ...
               'stinrf_export: %d labels but %d means', ...
               numel(fields.labels), numel(fields.mean));
    end
    if ~isempty(o.Matrices) && isfield(fields, 'mean')
        assert(numel(o.Matrices) == numel(fields.mean), ...
               'stinrf_export: %d matrices but %d means', ...
               numel(o.Matrices), numel(fields.mean));
    end

    d = fullfile(o.OutDir, name);
    if ~exist(d, 'dir'), mkdir(d); end

    for k = 1:numel(o.Matrices)
        writematrix(o.Matrices{k}, fullfile(d, sprintf('stim_%d.csv', k)));
    end
    if ~isempty(o.Curve)
        writematrix(o.Curve, fullfile(d, 'curve.csv'));
    end

    meta = fields;
    meta.name = string(name);
    meta.kind = string(kind);

    fid = fopen(fullfile(d, 'meta.json'), 'w');
    % Inf is not valid JSON; the plotter reads null as "show the full extent".
    fwrite(fid, strrep(jsonencode(meta, 'PrettyPrint', true), 'Inf', 'null'), 'char');
    fclose(fid);

    fprintf('exported %s\n', d);
end

function d = default_export_dir()
% DEFAULT_EXPORT_DIR  Where the export goes unless told otherwise.
    d = getenv('STINRF_EXPORT');
    if isempty(d)
        d = fullfile(inme_output(), 'export');
    end
end
