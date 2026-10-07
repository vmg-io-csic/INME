%% export_Opponency.m
%
% Writes the opponency figure data and a video of each stimulus. Runs in
% the workspace of run_Opponency_1D.m; add this line at the end of that
% script:
%
%     export_Opponency;
%
% No error bars (error_bars = false). make_Opponency_stimulus_1D takes no
% phase argument and draws nothing at random, so every trial builds the
% identical stimulus and the spread across trials is floating-point noise
% rather than variability. Setting n_trials to 1 in the run script would say
% the same thing and save the other 49 identical computations.

ex_n = 4;
ex_stims = cell(1, ex_n);  ex_labs = cell(1, ex_n);  ex_rlab = cell(1, ex_n);
ex_blab = cell(1, ex_n);
ex_dirs = zeros(1, ex_n);

for ex_i = 1:ex_n
    % Same pairing as the run loop: the first two conditions are the dark
    % bar, the last two the light one.
    if ex_i > 2, ex_inv = 1; else, ex_inv = 0; end
    ex_m = make_Opponency_stimulus_1D(N, n_frames, framerate, dov, ...
        sigma_t, speed_dps, bar_width, meanL, Direction(ex_i), ex_inv);
    ex_stims{ex_i} = meanL*(1 + 2*(ex_m./abs(max(max(ex_m))) - 0.5));

    if Direction(ex_i) >= 0, ex_word = 'Rightward'; else, ex_word = 'Leftward'; end
    % inver = 0 builds a bright bar on a black background and inver = 1
    % subtracts it from the mean, giving a dark bar on a light one. The
    % labels used to read the other way round, which the filmstrip hid and
    % the space-time panels do not: the panel called 'dark bar' showed a
    % white stripe. The responses are the same either way, so this is the
    % caption catching up with the stimulus.
    if ex_inv == 0, ex_bar = 'light bar'; else, ex_bar = 'dark bar'; end
    ex_labs{ex_i} = sprintf('%s, %s', ex_word, ex_bar);
    ex_rlab{ex_i} = {ex_word, ex_bar};
    % The bars get a short label, 'Light ->' and so on, which plot_figures
    % draws with an arrow: the panel above already says it in full.
    if Direction(ex_i) >= 0, ex_arr = '->'; else, ex_arr = '<-'; end
    ex_w = strtok(ex_bar);
    ex_blab{ex_i} = sprintf('%s %s', [upper(ex_w(1)) ex_w(2:end)], ex_arr);
    ex_dirs(ex_i) = sign(Direction(ex_i));
end

% Space-time, not a filmstrip: a single moving bar is one stripe per frame
% and the four conditions differ only in which way that stripe travels, which
% a filmstrip cannot show and a space-time plot shows as the slope of a line.
% The four panels go side by side in one row, with the bars beside them.
stinrf_export('Opponency', "panels", struct( ...
    'panel', "spacetime", 'space_win', Inf, 'frame_span', 0.65, 'ncols', 4, ...
    'dov', dov, 'N', N, 'framerate', framerate, 'n_frames', n_frames, ...
    'labels', {ex_labs}, 'bar_labels', {ex_blab}, 'rowlab', {ex_rlab}, 'dirs', ex_dirs, ...
    'mean', Mean_response(:).', 'sem', zeros(1, ex_n), 'error_bars', false), ...
    'Matrices', ex_stims);

% A video of each panel, from the same matrix. Grey levels are luminance:
% black is 0 and white twice the mean.
for ex_i = 1:ex_n
    stinrf_video(ex_stims{ex_i}, framerate, ...
        sprintf('Opponency_%d_%s', ex_i, ex_labs{ex_i}), 'Clim', [0 2*meanL]);
end

clear ex_n ex_stims ex_labs ex_rlab ex_dirs ex_i ex_inv ex_m ex_word ex_bar ...
      ex_blab ex_arr ex_w
