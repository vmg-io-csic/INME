%% export_MaskingNoMask.m
%
% Writes the signal-alone level that the masking figure draws as its
% reference line. Runs in the workspace of run_Masking_nomask_1D.m; add this
% line at the end of that script:
%
%     export_MaskingNoMask;
%
% This is not a figure of its own: plot_figures.py reads it while drawing
% Masking and does not draw it separately. Keeping it as an export rather
% than a constant in the masking script is what stops the reference drifting
% out of date when the model or the parameters change.

stinrf_export('MaskingNoMask', "reference", struct( ...
    'value', mean(Mean_response), 'signal_cpd', sfm_list(1)));

% A video of the signal alone, rebuilt with the run loop's call but at
% phase 0, with the same black and white as the masking videos so that the
% signal looks the same in all of them.
ex_m = make_Masking_stimulus_1D(N, n_frames, framerate, dov, sigma_t, 0, ...
    speed_dps, sfm_list, sfmask_cpd(1), meanL, masking, ...
    Contrast_signal, Contrast_mask, temporal_frequency_mask);
stinrf_video(ex_m, framerate, sprintf('MaskingNoMask_signal_%gcpd', sfm_list(1)), ...
    'Clim', meanL*(1 + (Contrast_signal + Contrast_mask)*[-1 1]));
clear ex_m
