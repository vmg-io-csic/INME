%% export_Contrast.m
%
% Writes the contrast response figure data and a video of the stimulus at
% every contrast. Runs in the workspace of run_Contrast_1D.m; add this line
% at the end of that script:
%
%     export_Contrast;
%
% The magnitude is exported, not the signed response: the figure is drawn on
% log axes, where a sign has nowhere to go.

stinrf_export('Contrast', "contrast", struct(), ...
    'Curve', [Contrasts(:), abs(Mean_response(:))]);

% A video at every contrast, the stimulus rebuilt exactly as the run loop
% builds it (phase fixed at 0). Grey levels are luminance, black 0 and white
% twice the mean, so the videos differ in contrast as the stimuli do.
for ex_i = 1:numel(Contrasts)
    ex_m = make_Contrast_stimulus_1D(N, n_frames, framerate, dov, ...
        sigma_t, 0, speed_dps, sfm_list, meanL);
    ex_m = meanL*(1 + Contrasts(ex_i)*2*(ex_m./abs(max(max(ex_m))) - 0.5));
    stinrf_video(ex_m, framerate, ...
        sprintf('Contrast_%02d_contrast_%.3f', ex_i, Contrasts(ex_i)), ...
        'Clim', [0 2*meanL]);
end
clear ex_i ex_m
