%% export_Gratings.m
%
% Writes the gratings figure. Runs in the workspace of run_Gratings.m; that
% script calls it on its last line.
%
% A grating has one component and no separate pattern direction, which is why
% its curve takes the component colour rather than the pattern one.

ex_icons = arrayfun(@(d) struct('comp', d, 'pattern', []), grating_dirs(:));

stinrf_export('Gratings', "directions", struct( ...
    'xlabel', "Stimulus direction (deg)", 'line_role', "comp1", ...
    'x', grating_dirs(:).', 'y', pred(:).', ...
    'mae', mean(abs(err)), 'worst', max(abs(err)), ...
    'icons', ex_icons));

clear ex_icons

% A video of the grating in every direction, rebuilt with the same call as
% the run loop. Grey levels are luminance: black is 0 and white twice the
% mean.
for ex_i = 1:nG
    ex_m = make_grating(N, n_frames, f_cpx, sp, grating_dirs(ex_i), phase, ...
                        contrast, meanL);
    stinrf_video(ex_m, framerate, ...
        sprintf('Gratings_%d_%gdeg', ex_i, grating_dirs(ex_i)), ...
        'Clim', meanL*(1 + contrast*[-1 1]));
end
clear ex_i ex_m
