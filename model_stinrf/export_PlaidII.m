%% export_PlaidII.m
%
% Writes the Type-II plaid figure. Runs in the workspace of run_PlaidII.m;
% that script calls it on its last line.
%
% Sorted by IOC direction so the axis runs left to right; the configurations
% are not stored in that order.

[ex_ioc, ex_ord] = sort(ioc);

ex_icons = arrayfun(@(i) struct( ...
    'comp', [plaidII_configs(ex_ord(i),3), plaidII_configs(ex_ord(i),4)], ...
    'pattern', ex_ioc(i)), (1:numel(ex_ioc)).');

stinrf_export('PlaidII', "directions", struct( ...
    'xlabel', "Stimulus direction (deg)", 'line_role', "pattern", ...
    'x', ex_ioc(:).', 'y', pred(ex_ord).', ...
    'mae', mean(abs(err)), 'worst', max(abs(err)), ...
    'icons', ex_icons));

clear ex_ioc ex_ord ex_icons

% A video of every configuration, in the order the run stores them, rebuilt
% with the same call as the run loop. The two components add, so the
% luminance spans meanL*(1 +/- 2c) and that range is black to white.
for ex_i = 1:nQ
    ex_m = make_plaidII(N, n_frames, f_cpx, sp1(ex_i), sp2(ex_i), ...
                        th1(ex_i), th2(ex_i), phase, phase, contrast, meanL);
    stinrf_video(ex_m, framerate, ...
        sprintf('PlaidII_%d_IOC_%gdeg', ex_i, round(ioc(ex_i), 1)), ...
        'Clim', meanL*(1 + 2*contrast*[-1 1]));
end
clear ex_i ex_m
