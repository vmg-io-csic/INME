%% run_all_experiments.m
%
% Runs every experiment in turn. Each one saves its own results file and,
% at the end, exports what its figure needs. There is no separate export
% step: an experiment that has run has exported.
%
% This is the slow path, for rebuilding everything from nothing. To redo a
% single figure, run its own script instead and redraw just that one:
%
%     matlab -batch "run_ReversePhi_1D"
%     python3 plot_figures.py <export folder> ReversePhi
%
% Outputs:
%   export/<Figure>/ for every figure, videos/<Figure>_<k>_<label>.mp4 for
%   every stimulus, and results/<Experiment>/ for the results files
%

clear; clc;

scripts = { ...
    'run_Masking_nomask_1D', ...   % before run_Masking_1D: it is the reference
    'run_Masking_1D', ...
    'run_Contrast_1D', ...
    'run_ReversePhi_1D', ...
    'run_Missing_1D', ...
    'run_SecondOrder_1D', ...
    'run_Opponency_1D', ...
    'run_Interaction_1D', ...
    'run_Gratings', ...
    'run_PlaidI', ...
    'run_PlaidII', ...
    'run_Uneven', ...
    'run_LightEdges'};

failed = {};
for k = 1:numel(scripts)
    fprintf('\n===== %s =====\n', scripts{k});
    try
        % Called through a function so that the script's own clear cannot
        % wipe the loop variables out from under this driver.
        stinrf_run_one(scripts{k});
    catch ME
        fprintf('  FAILED: %s\n', ME.message);
        failed{end+1} = scripts{k}; %#ok<SAGROW>
    end
end

fprintf('\n===== done, %d of %d ran =====\n', ...
        numel(scripts) - numel(failed), numel(scripts));
for k = 1:numel(failed), fprintf('  failed: %s\n', failed{k}); end
fprintf('\nThen draw everything with:\n  python3 plot_figures.py\n\n');
