%% export_Uneven.m
%
% Writes the uneven square wave figure and a video of each stimulus. Runs in the
% workspace of run_Uneven.m; that script calls it on its last line.
%
% The mean, the SD across phases (the whisker the figure draws) and the SEM
% are computed here, as in every other experiment, and written into
% meta.json. The per-phase responses go alongside as trials.csv
% so the summary can be checked against them, but nothing downstream
% recomputes it: MATLAB measures, Python draws.

ex_mu  = [mean(Rf_bar), mean(Rr_bar)];
ex_sd  = [std(Rf_bar), std(Rr_bar)];
ex_sem = ex_sd ./ sqrt(n_phases);

stinrf_export('Uneven', "frames_pair", struct( ...
    'title', "Uneven square wave", 'frame_span', 0.25, ...
    'dov', dov, 'N', N, 'framerate', framerate, ...
    'labels', {{'Forward', 'Reversed'}}, 'dirs', [1, -1], ...
    'n_trials', n_phases, 'mean', ex_mu, 'sem', ex_sem, 'sd', ex_sd), ...
    'Matrices', {S_stim, Sr_stim});

stinrf_export_bar('Uneven', [(1:n_phases)', Rf_bar(:), Rr_bar(:)]);

fprintf('  mean %+.4f / %+.4f, SD %.4f / %.4f, SEM %.4f / %.4f, n=%d\n', ...
        ex_mu(1), ex_mu(2), ex_sd(1), ex_sd(2), ex_sem(1), ex_sem(2), n_phases);

% A video of the forward stimulus and one of its time reverse, from the same
% matrices as the panels. Grey levels are luminance: black is 0 and white
% twice the mean. Built at 1000 fps, so real-time playback shows one frame
% in nine; 'Speed', 0.12 shows every one.
stinrf_video(S_stim,  framerate, 'Uneven_1_Forward',  'Clim', [0 2*meanL]);
stinrf_video(Sr_stim, framerate, 'Uneven_2_Reversed', 'Clim', [0 2*meanL]);

clear ex_mu ex_sem ex_sd
