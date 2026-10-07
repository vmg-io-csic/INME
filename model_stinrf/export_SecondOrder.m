%% export_SecondOrder.m
%
% Writes the second-order figure data and a video of each stimulus. Runs in
% the workspace of run_SecondOrder_1D.m; add this line at the end of that
% script:
%
%     export_SecondOrder;
%
% The whiskers are the SD across the n_trials runs.

ex_n = n_type;
ex_stims = cell(1, ex_n);  ex_labs = cell(1, ex_n);  ex_rlab = cell(1, ex_n);
ex_dirs = zeros(1, ex_n);

for ex_i = 1:ex_n
    % Both phases at 0, envelope and carrier, so the panel is reproducible.
    ex_stims{ex_i} = make_SecondOrder_stimulus_1D(N, n_frames, framerate, dov, ...
        sigma_t, 0, 0, speed_dps, sfm_envelope_cpd, sfm_carrier_cpd, ...
        meanL, Type(ex_i), Carrier_type);
    if Type(ex_i) < 0, ex_labs{ex_i} = 'Leftward'; else, ex_labs{ex_i} = 'Rightward'; end
    ex_rlab{ex_i} = {'Second order', ex_labs{ex_i}};
    ex_dirs(ex_i) = Type(ex_i);
end

% Space-time: what moves is the contrast envelope, and over five frames it
% travels too little to see. In a space-time plot it reads as a diagonal band.
stinrf_export('SecondOrder', "panels", struct( ...
    'panel', "spacetime", 'space_win', 1.25, ...
    'dov', dov, 'N', N, 'framerate', framerate, 'n_frames', n_frames, ...
    'labels', {ex_labs}, 'rowlab', {ex_rlab}, 'dirs', ex_dirs, ...
    'mean', Mean_response(:).', 'sem', SEM_response, 'sd', SD_response), ...
    'Matrices', ex_stims);

% A video of each panel, from the same matrix. Grey levels are luminance:
% black is 0 and white twice the mean.
for ex_i = 1:ex_n
    stinrf_video(ex_stims{ex_i}, framerate, ...
        sprintf('SecondOrder_%d_%s', ex_i, ex_labs{ex_i}), 'Clim', [0 2*meanL]);
end

clear ex_n ex_stims ex_labs ex_rlab ex_dirs ex_i

% The response of every run, kept beside meta.json so the mean and the SD
% can be checked against them. The figure does not draw them.
stinrf_export_bar('SecondOrder', [(1:n_trials)', All_trials.']);
