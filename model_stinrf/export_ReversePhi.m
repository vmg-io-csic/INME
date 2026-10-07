%% export_ReversePhi.m
%
% Writes the reverse-phi figure data and a video of each stimulus. This is a
% script, not a function, so it runs in the workspace of run_ReversePhi_1D.m
% and uses the parameters that experiment actually ran with. Add this line
% at the end of that script:
%
%     export_ReversePhi;
%
% The panels are drawn at phase 0 so the figure is reproducible; the bars
% still come from the randomised trials, with the SD across them as the
% whisker.

ex_n = n_type;
ex_stims = cell(1, ex_n);  ex_labs = cell(1, ex_n);  ex_rlab = cell(1, ex_n);
ex_dirs = zeros(1, ex_n);
if speed_dps >= 0, ex_word = 'Rightward'; else, ex_word = 'Leftward'; end

for ex_i = 1:ex_n
    ex_stims{ex_i} = make_ReversePhi_stimulus_1D(N, n_frames, framerate, dov, ...
        sigma_t, 0, speed_dps, meanL, Type(ex_i), Contrast);
    if Type(ex_i) > 0, ex_labs{ex_i} = 'Phi'; else, ex_labs{ex_i} = 'Reverse phi'; end
    ex_rlab{ex_i} = {ex_labs{ex_i}, ex_word};
    ex_dirs(ex_i) = sign(speed_dps);
end

% Space-time, not a filmstrip: the noise element is one pixel and the pattern
% steps four, so no single frame of reverse-phi differs from a phi frame and
% the correspondence between frames cannot be followed by eye. The diagonal
% streaks and the polarity alternation are obvious in a space-time plot.
stinrf_export('ReversePhi', "panels", struct( ...
    'panel', "spacetime", 'space_win', 2.0, ...
    'dov', dov, 'N', N, 'framerate', framerate, 'n_frames', n_frames, ...
    'labels', {ex_labs}, 'rowlab', {ex_rlab}, 'dirs', ex_dirs, ...
    'mean', Mean_response(:).', 'sem', SEM_response, 'sd', SD_response), ...
    'Matrices', ex_stims);

% A video of each panel, from the same matrix. Grey levels are luminance:
% black is 0 and white twice the mean. The polarity flips on every frame,
% so reverse phi only looks right on a screen that shows all 120 of them.
for ex_i = 1:ex_n
    stinrf_video(ex_stims{ex_i}, framerate, ...
        sprintf('ReversePhi_%d_%s', ex_i, ex_labs{ex_i}), 'Clim', [0 2*meanL]);
end

clear ex_n ex_stims ex_labs ex_rlab ex_dirs ex_word ex_i

% The response of every run, kept beside meta.json so the mean and the SD
% can be checked against them. The figure does not draw them.
stinrf_export_bar('ReversePhi', [(1:n_trials)', All_trials.']);
