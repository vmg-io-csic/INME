%% export_Masking.m
%
% Writes the masking figure data and a video of the signal with each mask.
% Runs in the workspace of run_Masking_1D.m; add this line at the end of
% that script:
%
%     export_Masking;
%
% The no-mask reference is not taken from the Response_no_mask constant in
% the run script: that value was measured under an older framerate and
% readout and is stale. It comes from the dedicated no-mask run instead, via
% its own export, so the two are always measured the same way.

stinrf_export('Masking', "masking", struct( ...
    'signal_cpd', sfm_list(1)), ...
    'Curve', [sfmask_cpd(:), Mean_response(:)]);

% A video of the signal with each mask, rebuilt with the run loop's call but
% at phase 0. The mask jitters at random (make_Masking_stimulus_1D draws its
% direction afresh on every frame), so each video is one draw of it. Black
% and white are the furthest signal and mask can reach together.
ex_clim = meanL*(1 + (Contrast_signal + Contrast_mask)*[-1 1]);
for ex_i = 1:n_masking
    ex_m = make_Masking_stimulus_1D(N, n_frames, framerate, dov, sigma_t, 0, ...
        speed_dps, sfm_list, sfmask_cpd(ex_i), meanL, masking, ...
        Contrast_signal, Contrast_mask, temporal_frequency_mask);
    stinrf_video(ex_m, framerate, ...
        sprintf('Masking_%02d_mask_%gcpd', ex_i, sfmask_cpd(ex_i)), 'Clim', ex_clim);
end
clear ex_i ex_m ex_clim
