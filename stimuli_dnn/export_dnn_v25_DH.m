%% export_dnn_v25_DH.m
%
% Exports the Derrington & Henning / interaction-among-scales stimuli for
% the five optical-flow networks. There was no DNN export for DH before.
%
% TWO CHANGES, both sampling only:
%   px_per_deg 100 -> 180.  DH is the one family where this is forced: with
%       ppd = 100 the 3 cpd moving component tops out at 8.3 px/frame and
%       cannot reach the band whatever the framerate. At 180 both
%       components land inside it, and dov becomes 4.27 deg -- essentially
%       the dov = 4.0 that stINRF3_config declares for DH (the old export
%       was at 7.68 deg).
%   framerate 1000 -> 60.   Gives 0.20 cycles/frame.
%
% The stimulus is untouched: 1 and 3 cpd, 12 Hz, the four sigma_t.
% Duration stays at stINRF3's 500 ms so the Gaussian envelope occupies the
% same fraction of the clip in both models -> 30 frames.
%
%   1 cpd (S1, CS)      12 deg/s  -> 36 px/frame
%   3 cpd (S3, CSr, CM)  4 deg/s  -> 12 px/frame
%   both at 0.20 cycles/frame
%
% Uses make_DH_stimulus_1D_v25 (bug B2 fixed: the drift now depends on
% `framerate`, and the misnamed `speed_dps` is `tf_hz`).
%
% CAVEAT worth reporting: at sigma_t = 10 ms the stimulus is visible for
% 23.6 ms, i.e. 2 frames at 60 fps. That is enough for RAFT, PWC-Net and
% LiteFlowNet2 (they consume pairs) but ME-Attention and Dual receive a
% window in which only 2 of 11 / 15 frames carry signal. No sampling rate
% fixes this: 15 frames inside 23.6 ms would need 636 fps and a ppd so high
% that less than one cycle of the 1 cpd component fits on the canvas.

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

meanL = 100;  phase_mov = 0;  phase_stat = 0;

for g = 1:numel(groups)
    tag = groups(g).tag;  N_cols = groups(g).N_cols;  N_rows = groups(g).N_rows;
    c = dnn_v25_config('DH', N_cols);
    out_root = fullfile('stimuli_by_model', tag, 'v25');
    fprintf('=== DH %s (%dx%d) | ppd=%g fr=%g dov=%.2f deg | %d frames ===\n', ...
            tag, N_cols, N_rows, c.px_per_deg, c.framerate, c.dov, c.n_frames);
    count = 0; skipped = 0;

    for ci = 1:numel(c.cond_names)
        for si = 1:numel(c.sigma_t_ms)
            sigma_t = c.sigma_t_ms(si);
            name = sprintf('DH_%s_sig%05.1f', c.cond_names{ci}, sigma_t);
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file')
                skipped = skipped + 1; count = count + 1; continue;
            end

            mI = make_DH_stimulus_1D_v25(N_cols, c.n_frames, c.framerate, ...
                    c.dov, ci, sigma_t, phase_mov, phase_stat, c.tf_hz, ...
                    c.sfm, c.sfs, meanL);

            sfm = c.sfm(ci);
            v_dps = c.tf_hz / sfm;
            px_f  = v_dps * c.px_per_deg / c.framerate;
            % how much of the window actually carries signal
            fwhm_ms = 2.355 * sigma_t;
            fr_signal = fwhm_ms / 1000 * c.framerate;

            extra = struct('experiment','DH_v25', ...
                'condition', c.cond_names{ci}, 'cond_idx', ci, ...
                'sfm_cpd', sfm, 'sfs_cpd', c.sfs(ci), ...
                'tf_hz', c.tf_hz, 'speed_dps', v_dps, ...
                'sigma_t_ms', sigma_t, ...
                'px_per_deg', c.px_per_deg, 'framerate', c.framerate, ...
                'dov', c.dov, 'N_cols', N_cols, 'N_rows', N_rows, ...
                'meanL', meanL, ...
                'speed_px_per_frame', px_f, ...
                'cycles_per_frame', c.tf_hz / c.framerate, ...
                'frames_within_fwhm', fr_signal, ...
                'gt_direction_deg', 0);       % rightward by construction

            dnn_v25_save(edir, dnn_v25_broadcast(mI, N_rows), extra);
            count = count + 1;
            if si == 1
                fprintf('  %-4s sfm=%d cpd  v=%5.2f deg/s  %5.1f px/frame\n', ...
                        c.cond_names{ci}, sfm, v_dps, px_f);
            end
        end
    end
    fprintf('  %d stimuli (%d skipped)\n\n', count, skipped);
end
fprintf('Done.\n');
