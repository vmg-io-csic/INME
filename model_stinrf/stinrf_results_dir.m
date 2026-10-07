function d = stinrf_results_dir(script_name)
% STINRF_RESULTS_DIR  Where an experiment saves its results file.
%
%   One folder per experiment, under results/, rather than every run of every
%   experiment in a single directory. With everything in one place a stale
%   file from weeks earlier can be picked up in place of the one just
%   written, which is easy to miss because the figure still draws.
%
%   Results files are the record of a run. The figures do not read them: they
%   read the export.
%
% Inputs
%   script_name (char) usually mfilename, for example 'run_ReversePhi_1D'
%
% Outputs
%   d (char) the folder, created if it does not exist

    name = char(script_name);
    name = regexprep(name, '^(run|gen)_', '');
    name = regexprep(name, '_1D$|_2D$', '');
    d = fullfile(inme_output(), 'results', name);
    if ~exist(d, 'dir'), mkdir(d); end
end
