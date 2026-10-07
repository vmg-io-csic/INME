%% export_Interaction.m
%
% Writes the cross-scale interaction figure and a video of each stimulus.
% Runs in the workspace of run_Interaction_1D.m; add this line at the end
% of that script:
%
%     export_Interaction;
%
% The experiment measures every condition at every envelope width, and the
% results file keeps all of it. The figure shows one condition at three
% widths, because what it has to make visible is that the reported direction
% reverses as the envelope widens, and a panel per width says that where a
% curve does not. Everything else stays in the .mat for the record.

ex_cond  = 4;                      % CSr: 3 cpd moving on a 1 cpd static
ex_ends  = [10, 200];              % short and long, well either side

ex_ic = find(cond_subset == ex_cond, 1);
if isempty(ex_ic)
    error(['export_Interaction: condition %d was not run; cond_subset is ' ...
           'currently [%s]'], ex_cond, num2str(cond_subset));
end

ex_mu_all  = mean(val, 3);
ex_sd_all  = std(val, 0, 3);
ex_sem_all = ex_sd_all ./ sqrt(n_trials);

% The figure plots the share of trials the model got right rather than the
% mean response. The responses here are around a hundredth of a unit, so a
% bar chart in arbitrary units says nothing until the reader has worked out
% which sign means which way; the share of trials whose sign matches the
% direction the stimulus actually drifted in says it outright. Same
% convention as pct_correct in the run script.
ex_pc_all  = mean(val * direction > 0, 3);

% The middle panel is the width at which the model is closest to guessing,
% read off the sweep rather than written in by hand. It used to be the
% constant 89.2 ms, which was where an earlier run had put the crossing;
% once the sweep moved, the panel labelled "at the crossing" was at 64 per
% cent correct and the label was no longer true. Taking the measured
% minimum keeps the panel honest whatever the sweep does next.
[~, ex_im] = min(abs(ex_pc_all(ex_ic, :) - 0.5));
ex_sigma = [ex_ends(1), sigma_t_list(ex_im), ex_ends(2)];
fprintf('  closest to chance: %g ms at %.0f%% correct\n', ...
        sigma_t_list(ex_im), 100*ex_pc_all(ex_ic, ex_im));

ex_n = numel(ex_sigma);
ex_stims = cell(1, ex_n);  ex_mu = zeros(1, ex_n);  ex_sem = zeros(1, ex_n);
ex_caps = cell(1, ex_n);   ex_labs = cell(1, ex_n);  ex_pc = zeros(1, ex_n);
ex_sd = zeros(1, ex_n);

for ex_i = 1:ex_n
    [ex_d, ex_is] = min(abs(sigma_t_list - ex_sigma(ex_i)));
    if ex_d > 1e-6
        warning(['export_Interaction: no envelope width at %g ms, using ' ...
                 '%g ms instead'], ex_sigma(ex_i), sigma_t_list(ex_is));
    end
    % Both phases at 0, so the panel is reproducible; the bars still come
    % from the randomised trials the experiment ran.
    ex_stims{ex_i} = make_Interaction_stimulus_1D(N, n_frames, framerate, dov, ...
        cond_subset(ex_ic), sigma_t_list(ex_is), 0, 0, ...
        v_High, v_Low, sfm_High_cpd, sfm_Low_cpd, meanL, direction);
    ex_mu(ex_i)  = ex_mu_all(ex_ic, ex_is);
    ex_sem(ex_i) = ex_sem_all(ex_ic, ex_is);
    ex_sd(ex_i)  = ex_sd_all(ex_ic, ex_is);
    ex_pc(ex_i)  = ex_pc_all(ex_ic, ex_is);
    ex_caps{ex_i} = sprintf('sigma_t = %g ms', sigma_t_list(ex_is));
    ex_labs{ex_i} = sprintf('%g ms', sigma_t_list(ex_is));
end

% The stimulus drifts the same way in every panel, so the direction is stated
% once above the row rather than carried in the colour of each frame.
ex_bar_dir = sign(ex_mu);  ex_bar_dir(ex_bar_dir == 0) = 1;
if direction >= 0, ex_word = "Rightward"; else, ex_word = "Leftward"; end

% Binomial standard error of the share correct, kept in meta.json for the
% record. The bars are drawn without whiskers: the spread of a share correct
% is fixed by the share itself (it is largest at 50 per cent whatever the
% model does), so it would say nothing the bar does not.
ex_pc_sem = sqrt(ex_pc .* (1 - ex_pc) ./ n_trials);
% What the figure shows instead is the band of shares a model answering at
% chance would produce: every k out of n_trials correct that a two-sided
% exact binomial test against p = 0.5 does not reject at alpha = 0.05. With
% 25 trials it is 8 to 17 correct, 32 to 68 per cent. A bar inside the band
% cannot be told from guessing. Summed by hand from the binomial pmf so it
% needs no toolbox.
ex_k   = 0:n_trials;
ex_pmf = exp(gammaln(n_trials+1) - gammaln(ex_k+1) - gammaln(n_trials-ex_k+1) ...
             - n_trials*log(2));
ex_cdf = cumsum(ex_pmf);                         % P(X <= k)
ex_sf  = fliplr(cumsum(fliplr(ex_pmf)));         % P(X >= k)
ex_p2  = min(1, 2*min(ex_cdf, ex_sf));           % two-sided p of each k
ex_kin = ex_k(ex_p2 >= 0.05);
ex_band = [min(ex_kin), max(ex_kin)] / n_trials;
fprintf('  chance band %.0f-%.0f%% correct (n = %d)\n', 100*ex_band, n_trials);

stinrf_export('Interaction', "st_row", struct( ...
    'dov', dov, 'N', N, 'framerate', framerate, 'n_frames', n_frames, ...
    'space_win', 1.0, ...
    'n_trials', n_trials, 'condition', string(cond_names{ex_ic}), ...
    'captions', {ex_caps}, 'labels', {ex_labs}, ...
    'panel_dirs', direction*ones(1, ex_n), 'dirs', ex_bar_dir, ...
    'mean', ex_mu, 'sem', ex_sem, 'sd', ex_sd, ...
    'correct', ex_pc, 'correct_sem', ex_pc_sem, ...
    'chance', 0.5, 'chance_band', ex_band, ...
    'direction_word', ex_word), ...
    'Matrices', ex_stims);

fprintf('  figure condition %s at %s ms\n', cond_names{ex_ic}, ...
        strjoin(ex_labs, ', '));

fprintf('  correct %s\n', strjoin(compose('%.0f%%', 100*ex_pc), ', '));

% A video of each panel, from the same matrix. Grey levels are luminance:
% black is 0 and white twice the mean. The stimulus is built at 1000 fps,
% so real-time playback shows one frame in nine; the 10 ms envelope then
% lasts only a few frames, and 'Speed', 0.12 shows every one of them.
for ex_i = 1:ex_n
    stinrf_video(ex_stims{ex_i}, framerate, ...
        sprintf('Interaction_%d_sigma_t_%s', ex_i, ex_labs{ex_i}), ...
        'Clim', [0 2*meanL]);
end

clear ex_cond ex_ends ex_sigma ex_im ex_ic ex_mu_all ex_sem_all ex_pc_all ex_n ex_stims ...
      ex_mu ex_sem ex_pc ex_pc_sem ex_caps ex_labs ex_i ex_d ex_is ...
      ex_sd_all ex_sd ex_k ex_pmf ex_cdf ex_sf ex_p2 ex_kin ex_band ...
      ex_bar_dir ex_word
