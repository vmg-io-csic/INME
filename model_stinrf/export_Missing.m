%% export_Missing.m
%
% Writes the missing-fundamental figure data and a video of each stimulus.
% Runs in the workspace of run_Missing_1D.m; add this line at the end of
% that script:
%
%     export_Missing;
%
% No error bars: this is a single run at a fixed phase (the response changes
% by less than 0.2% across phases, see run_Missing_1D.m), so there is no
% spread across runs to show. error_bars = false tells the figure so.

ex_n = n_type;
ex_stims = cell(1, ex_n);  ex_labs = cell(1, ex_n);  ex_rlab = cell(1, ex_n);
ex_dirs = zeros(1, ex_n);
if speed_dps >= 0, ex_word = 'Rightward'; else, ex_word = 'Leftward'; end

for ex_i = 1:ex_n
    ex_stims{ex_i} = make_Missing_stimulus_1D(N, n_frames, framerate, dov, ...
        sigma_t, 0, speed_dps, sfm_list, meanL, Type(ex_i));
    if Type(ex_i) > 0
        ex_labs{ex_i} = 'Square wave';
    else
        ex_labs{ex_i} = 'Missing fundamental';
    end
    ex_rlab{ex_i} = {ex_labs{ex_i}, ex_word};
    ex_dirs(ex_i) = sign(speed_dps);
end

% Space-time, not a filmstrip. Five frames show five instants and leave the
% shift between them to be inferred; in a space-time plot the square wave and
% the missing fundamental both read as diagonal bands and the reader can see
% that the two patterns drift the same way while the model reports opposite
% directions, which is the whole result.
stinrf_export('Missing', "panels", struct( ...
    'panel', "spacetime", 'space_win', 1.25, 'frame_span', 0.20, ...
    'dov', dov, 'N', N, 'framerate', framerate, 'n_frames', n_frames, ...
    'labels', {ex_labs}, 'rowlab', {ex_rlab}, 'dirs', ex_dirs, ...
    'mean', Mean_response(:).', 'sem', SEM_response, 'error_bars', false), ...
    'Matrices', ex_stims);

% A video of each panel, from the same matrix. Both patterns run past
% [0, 2*meanL] (to about -13 and 213), so black and white are the lowest and
% highest values of the two, as in the figure. The range is symmetric about
% the mean, which stays mid grey.
ex_clim = [min(cellfun(@(m) min(m(:)), ex_stims)), ...
           max(cellfun(@(m) max(m(:)), ex_stims))];
for ex_i = 1:ex_n
    stinrf_video(ex_stims{ex_i}, framerate, ...
        sprintf('Missing_%d_%s', ex_i, ex_labs{ex_i}), 'Clim', ex_clim);
end

clear ex_n ex_stims ex_labs ex_rlab ex_dirs ex_word ex_i ex_clim
