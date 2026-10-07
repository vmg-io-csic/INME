%% export_dnn_v25_gratingmix.m
%
% Grating plus co-moving texture, with the TEXTURE held fixed and the
% GRATING contrast swept down to zero.
%
% WHY
% ---
% export_dnn_v25_gratingtex.m sweeps the texture upwards against a fixed
% grating, and on that ladder RAFT is fixed by the very first step while
% ME-Attention does not respond at all and Dual only improves gradually.
%
% The hypothesis for that split is about which component owns the contrast
% energy. RAFT matches: any distinctive structure gives a correspondence,
% so a weak texture is enough. ME-Attention and Dual are motion-energy
% models built from banks of oriented spatiotemporal filters, and their
% population response is dominated by whatever carries the most contrast.
% On the gratingTex ladder the grating still holds about 70% of the energy
% at the top step (RMS 0.283 against 0.12), so those two keep seeing what
% is essentially a grating.
%
% This file tests that directly by removing the grating's energy instead
% of adding to the texture's. The texture is fixed at 0.12 RMS and the
% grating contrast goes 0.2, 0.1, 0.0, which puts the grating at 58%, 26%
% and 0% of the total energy. c_grating = 0.4 is not repeated here: it is
% already gratingTex_ct120, the top of the other ladder.
%
% The prediction is explicit, so the experiment can fail: ME-Attention and
% Dual should improve sharply as the grating's share falls, and RAFT should
% not change, since it already succeeds at 70%.
%
% c_grating = 0 is also the positive control the design otherwise lacks --
% pure translating texture, no grating anywhere. A model that fails there
% is not failing because of the grating.
%
% SAMPLING
% --------
% dnn_v25_config('texture'): 85.333 px/deg, 85.33 fps, 36 frames,
% 16.00 px/frame, 0.5 cpd, directions 0:15:90, texture band 40-60 px.
%
% Resume-safe. Set SWEEP_ONLY_TAG=<tag> for one canvas group.

clear; clc;

groups_all = struct( ...
    'tag',    {'raft',  'pwc_lfn2', 'meattn_dual'}, ...
    'N_cols', {960,     1024,       768}, ...
    'N_rows', {520,     436,        768});

only_tag = getenv('SWEEP_ONLY_TAG');
if ~isempty(only_tag)
    idx = find(strcmp({groups_all.tag}, only_tag));
    if isempty(idx), error('SWEEP_ONLY_TAG=%s unknown', only_tag); end
    groups = groups_all(idx);
else
    groups = groups_all;
end

meanL      = 100;
c_texture  = 0.12;
c_grat_list = [0.2 0.1 0.0];
period_lo  = 40;
period_hi  = 60;
seed_tex   = 2;
phase0     = 0;

for g = 1:numel(groups)
    tag = groups(g).tag;  N_cols = groups(g).N_cols;  N_rows = groups(g).N_rows;
    out_root = fullfile('stimuli_by_model', tag, 'v25');
    c = dnn_v25_config('texture', N_cols);
    f_cpx = c.f_spatial / c.px_per_deg;
    speed = c.speed_Hz / c.framerate;
    v_px  = speed / f_cpx;
    fprintf('=== %s (%dx%d)  v = %.3f px/frame ===\n', tag, N_cols, N_rows, v_px);

    count = 0; skipped = 0;
    for cg = c_grat_list
        for d = c.directions_deg
            name = sprintf('gratingMix_dir%03d_cg%03d', d, round(cg*1000));
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file')
                skipped = skipped + 1;  count = count + 1;  continue;
            end
            [TCB, info] = make_textured_grating(N_cols, c.n_frames, f_cpx, ...
                    speed, d, phase0, cg, c_texture, meanL, ...
                    period_lo, period_hi, seed_tex);
            if info.clipped
                error('gratingMix %s clipped', name);
            end
            rms_g = cg / sqrt(2);
            extra = struct('experiment','gratingMix_v25', ...
                'direction_deg', d, 'gt_direction_deg', d, ...
                'grating_contrast', cg, 'texture_rms_contrast', c_texture, ...
                'grating_energy_fraction', rms_g^2/(rms_g^2 + c_texture^2), ...
                'texture_period_lo_px', period_lo, ...
                'texture_period_hi_px', period_hi, 'seed_texture', seed_tex, ...
                'f_spatial_cpd', c.f_spatial, 'speed_Hz', c.speed_Hz, ...
                'px_per_deg', c.px_per_deg, 'framerate', c.framerate, ...
                'dov', c.dov, 'N_cols', N_cols, 'N_rows', N_rows, ...
                'meanL', meanL, 'speed_px_per_frame', v_px, ...
                'cycles_per_frame', c.speed_Hz / c.framerate, ...
                'max_excursion', info.max_excursion);
            dnn_v25_save(edir, dnn_v25_crop(TCB, N_rows), extra);
            count = count + 1;
        end
    end
    fprintf('  gratingMix: %d (%d ya existian)\n', count, skipped);
end
fprintf('Done.\n');
