%% export_dnn_v25_phenomena.m
%
% Exports Contrast, Masking, Missing Fundamental, Opponency, Reverse Phi,
% Second Order and Time Reversal for the five optical-flow networks,
% resampled into the 10-40 px/frame band.
%
% WHAT CHANGES, AND WHY
% ---------------------
% Sampling only, except where noted:
%
%   Contrast / Masking     ppd 100 unchanged. ONE framerate for the whole
%   / SecondOrder          experiment, so px/frame does not co-vary with sf
%                          (otherwise a difference across the sweep could
%                          not be attributed to either axis). The finest sf
%                          sets it:
%                            Contrast, Masking  82 fps -> 11.4 px/f, all sf
%                            SecondOrder        47 fps -> 20.0 px/f, all sf
%                          sf = 4 cpd is DROPPED: its ceiling is 6.25
%                          px/frame, it can never reach the band at ppd 100.
%
%   SecondOrder            carrier grain 1 px -> 16 px. A 1 px binary
%                          carrier does not survive the networks' x8
%                          downsampling, so the envelope would modulate a
%                          texture the network cannot see.
%
%   Missing                speed 9.38 -> 5.0 deg/s and framerate 25 fps
%                          (= one frame per 40 ms step). The stimulus moves
%                          in steps, so at any other rate the 2-frame
%                          models get pairs with ZERO displacement. sf is
%                          capped at 1 cpd: the harmonic series runs to
%                          39*sf and must stay under the pixel Nyquist
%                          (ppd/2 = 50 cpd), so the existing sf = 2 and 4
%                          conditions were already aliased.
%
%   Opponency              ppd 100 -> 200 so the narrowest bar (0.125 deg)
%                          is 25 px and can take a 12 px displacement.
%                          A single ppd for all four widths keeps them
%                          comparable. dov = 3.84 deg, close to stINRF3's
%                          4.0. Uses the framerate-correct generator
%                          (bug B1).
%
%   ReversePhi             noise grain 1 px -> 50 px (>= 2x displacement);
%                          white noise is uncorrelated beyond 1 px so a
%                          20 px jump has nothing to match. Generated
%                          directly at 11.7 fps -- it cannot be decimated,
%                          the effect lives in the frame-to-frame contrast
%                          inversion.
%
%   TimeReversal           velocities and periods of the paper untouched.
%                          LightEdges 40 fps, UnevenSquareWave 5 fps, both
%                          at 20 px/frame. The second one needs a 7.2 s
%                          clip to reach 36 frames: it is simply watched
%                          for longer, at 1.17 deg/s that is 3 periods.
%                          Both use make_TimeReversal_stimulus_1D unchanged.
%                          NOTE: LightEdges is deliberately NOT a rigid
%                          translation -- it is Wu's S(x,t) = sign(vt - x)
%                          on 0<=x<lambda, 0<=t<T with lambda = v*T,
%                          periodic in both. Its light region grows to fill
%                          the period and resets, so the frame-mean
%                          luminance sweeps the full range. That asymmetry
%                          is the point: a rigidly translating square wave
%                          possesses Gamma, chi.Theta and chi.Gamma.Theta
%                          and is therefore symmetric under time reversal
%                          (Wu, Appendix 2), which would erase the effect.
%
% All clips are 36 frames, so Dual (15-frame windows) gets 22 and
% ME-Attention (11) gets 26.
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

meanL = 100;  phase_mov = 0;  phase_carrier = 0;

for g = 1:numel(groups)
    tag = groups(g).tag;  N_cols = groups(g).N_cols;  N_rows = groups(g).N_rows;
    out_root = fullfile('stimuli_by_model', tag, 'v25');
    fprintf('=== %s (%dx%d) ===\n', tag, N_cols, N_rows);
    count = 0; skipped = 0;

    % helper closure for the common metadata
    base = @(ppd, fr, v) struct('px_per_deg',ppd,'framerate',fr, ...
        'dov',N_cols/ppd,'N_cols',N_cols,'N_rows',N_rows,'meanL',meanL, ...
        'speed_dps',v,'speed_px_per_frame', v*ppd/fr);

    %% ================= CONTRAST =================
    c = dnn_v25_config('Contrast', N_cols);
    for si = 1:numel(c.sf_list)
        sf = c.sf_list(si);  fr = c.framerate;
        dov = N_cols / c.px_per_deg;
        for ci = 1:numel(c.contrast_list)
            con = c.contrast_list(ci);
            name = sprintf('Contrast_sf%.3f_c%.3f', sf, con);
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
            mI = make_Contrast_stimulus_1D(N_cols, c.n_frames, fr, dov, ...
                    c.sigma_t_ms, phase_mov, c.speed_dps, sf, meanL);
            mI = meanL * (1 + con * 2 * (mI ./ abs(max(mI(:))) - 0.5));
            e = base(c.px_per_deg, fr, c.speed_dps);
            e.experiment='Contrast_v25'; e.sfm_cpd=sf; e.contrast_level=con;
            e.sigma_t_ms=c.sigma_t_ms; e.cycles_per_frame=sf*c.speed_dps/fr;
            e.gt_direction_deg=0;
            dnn_v25_save(edir, dnn_v25_broadcast(mI, N_rows), e);
            count = count + 1;
        end
    end
    fprintf('  Contrast: %d\n', count);

    %% ================= MASKING =================
    c = dnn_v25_config('Masking', N_cols);
    for si = 1:numel(c.sf_list)
        sf = c.sf_list(si);  fr = c.framerate;
        dov = N_cols / c.px_per_deg;
        % the mask drifts at 10 Hz; framerate >= 4*10 keeps it unaliased
        if fr < 4*c.tf_mask_hz
            warning('Masking sf=%.2f: fr=%.1f < 4*%.0f Hz mask', sf, fr, c.tf_mask_hz);
        end
        % vi = 0 is the unmasked control, vi >= 1 the three mask ratios
        for vi = 0:numel(c.mask_ratios)
            if vi == 0
                msk = 0; ratio = 0; sfmask = sf;
                name = sprintf('Masking_sf%.3f_nomask', sf);
            else
                msk = 1; ratio = c.mask_ratios(vi); sfmask = sf * ratio;
                name = sprintf('Masking_sf%.3f_maskr%.2f', sf, ratio);
            end
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
            mI = make_Masking_stimulus_1D(N_cols, c.n_frames, fr, dov, ...
                    c.sigma_t_ms, phase_mov, c.speed_dps, sf, sfmask, meanL, ...
                    msk, c.contrast_signal, c.contrast_mask, c.tf_mask_hz);
            e = base(c.px_per_deg, fr, c.speed_dps);
            e.experiment='Masking_v25'; e.masking=msk; e.sfm_cpd=sf;
            e.sfmask_cpd=sfmask; e.mask_ratio=ratio; e.sigma_t_ms=c.sigma_t_ms;
            e.cycles_per_frame=sf*c.speed_dps/fr;
            e.mask_cycles_per_frame=c.tf_mask_hz/fr; e.gt_direction_deg=0;
            dnn_v25_save(edir, dnn_v25_broadcast(mI, N_rows), e);
            count = count + 1;
        end
    end
    fprintf('  +Masking: %d\n', count);

    %% ================= MISSING FUNDAMENTAL =================
    c = dnn_v25_config('Missing', N_cols);
    dov = N_cols / c.px_per_deg;
    for sf = c.sf_list
        if 39*sf >= c.px_per_deg/2
            error('Missing sf=%.2f: harmonic 39*sf=%.1f cpd exceeds Nyquist %.1f', ...
                  sf, 39*sf, c.px_per_deg/2);
        end
        for ty = c.types
            name = sprintf('Missing_sf%.3f_type%+d', sf, ty);
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
            mI = make_Missing_stimulus_1D(N_cols, c.n_frames, c.framerate, ...
                    dov, c.sigma_t_ms, phase_mov, c.speed_dps, sf, meanL, ty);
            e = base(c.px_per_deg, c.framerate, c.speed_dps);
            e.experiment='Missing_v25'; e.type=ty; e.sfm_cpd=sf;
            e.sigma_t_ms=c.sigma_t_ms;
            e.step_ms=40; e.px_per_step=c.speed_dps*0.040*c.px_per_deg;
            e.gt_direction_deg = (ty>0)*0 + (ty<0)*180;
            dnn_v25_save(edir, dnn_v25_broadcast(mI, N_rows), e);
            count = count + 1;
        end
    end
    fprintf('  +Missing: %d\n', count);

    %% ================= OPPONENCY =================
    c = dnn_v25_config('Opponency', N_cols);
    dov = N_cols / c.px_per_deg;
    for bw = c.bar_widths
        disp_px = c.speed_dps * c.px_per_deg / c.framerate;
        if disp_px > bw * c.px_per_deg / 2
            warning('Opponency bw=%.3f deg: %.1f px/frame exceeds half the bar (%.1f px)', ...
                    bw, disp_px, bw*c.px_per_deg/2);
        end
        % Both contrast polarities, as in run_Opponency_1D.m:
        % inver = 1 barra oscura sobre fondo claro, inver = 0 al reves.
        for iv = c.inver
            for dr = c.directions
                name = sprintf('Opponency_bw%.3f_dir%+d_inv%d', bw, dr, iv);
                edir = fullfile(out_root, name);
                if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
                mI = make_Opponency_stimulus_1D_v25(N_cols, c.n_frames, ...
                        c.framerate, dov, c.sigma_t_ms, c.speed_dps, bw, meanL, ...
                        dr, iv);
                mI = meanL * (1 + 2 * (mI ./ abs(max(mI(:))) - 0.5));
                e = base(c.px_per_deg, c.framerate, c.speed_dps);
                e.experiment='Opponency_v25'; e.direction=dr; e.bar_width=bw;
                e.bar_width_px=bw*c.px_per_deg; e.inver=iv;
                e.polarity_desc = polarity_label(iv);
                e.sigma_t_ms=c.sigma_t_ms;
                e.gt_direction_deg = (dr>0)*0 + (dr<0)*180;
                dnn_v25_save(edir, dnn_v25_broadcast(mI, N_rows), e);
                count = count + 1;
            end
        end
    end
    fprintf('  +Opponency: %d\n', count);

    %% ================= REVERSE PHI =================
    c = dnn_v25_config('ReversePhi', N_cols);
    dov = N_cols / c.px_per_deg;
    for con = c.contrast_list
        for ty = c.types
            name = sprintf('ReversePhi_c%.3f_type%+d', con, ty);
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
            mI = make_ReversePhi_stimulus_1D_grain(N_cols, c.n_frames, ...
                    c.framerate, dov, c.sigma_t_ms, phase_mov, c.speed_dps, ...
                    meanL, ty, con, c.grain_px);
            e = base(c.px_per_deg, c.framerate, c.speed_dps);
            e.experiment='ReversePhi_v25'; e.type=ty; e.contrast=con;
            e.grain_px=c.grain_px; e.sigma_t_ms=c.sigma_t_ms;
            % type -1 inverts contrast at every step: direction reverses
            e.gt_direction_deg = (ty>0)*0 + (ty<0)*180;
            dnn_v25_save(edir, dnn_v25_broadcast(mI, N_rows), e);
            count = count + 1;
        end
    end
    fprintf('  +ReversePhi: %d\n', count);

    %% ================= SECOND ORDER =================
    c = dnn_v25_config('SecondOrder', N_cols);
    dov = N_cols / c.px_per_deg;
    for si = 1:numel(c.sf_list)
        sf = c.sf_list(si);  fr = c.framerate;
        for ty = [-1, 1]
            name = sprintf('SecondOrder_sf%.3f_type%+d', sf, ty);
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
            mI = make_SecondOrder_stimulus_1D_grain(N_cols, c.n_frames, fr, ...
                    dov, c.sigma_t_ms, phase_mov, phase_carrier, c.speed_dps, ...
                    sf, 8, meanL, ty, 1, c.carrier_grain_px);
            e = base(c.px_per_deg, fr, c.speed_dps);
            e.experiment='SecondOrder_v25'; e.type=ty; e.sfm_envelope_cpd=sf;
            e.carrier_grain_px=c.carrier_grain_px; e.carrier_type=1;
            e.sigma_t_ms=c.sigma_t_ms; e.cycles_per_frame=sf*c.speed_dps/fr;
            e.gt_direction_deg = (ty>0)*0 + (ty<0)*180;
            dnn_v25_save(edir, dnn_v25_broadcast(mI, N_rows), e);
            count = count + 1;
        end
    end
    fprintf('  +SecondOrder: %d\n', count);

    %% ================= TIME REVERSAL =================
    % Excluded by default: the sampling density of the stimulus' own cycle
    % (T = lambda/v = 125 ms) is still open -- at ppd 85.33 / 40 fps only 5
    % frames fall in each cycle, so 1 in 5 pairs straddles the reset.
    % Run with  DNN_V25_TIMEREVERSAL=1  once that is decided.
    if ~strcmp(getenv('DNN_V25_TIMEREVERSAL'), '1')
        fprintf('  TimeReversal: SKIPPED (set DNN_V25_TIMEREVERSAL=1 to include)\n');
        fprintf('  %d stimuli total (%d skipped)\n\n', count, skipped);
        continue
    end
    c = dnn_v25_config('TimeReversal', N_cols);
    dov = N_cols / c.px_per_deg;
    for k = 1:numel(c.stim)
        s = c.stim(k);
        mI0 = make_TimeReversal_stimulus_1D(s.name, N_cols, c.n_frames, ...
                s.framerate, dov, 0, 1, meanL, s.v_dps, s.period_deg);
        for rev = [0 1]
            mI = mI0;
            if rev, mI = flipud(mI); end   % time reversal, NOT sign flip
            if rev, name = sprintf('TR_%s_reversed', s.name);
            else,   name = sprintf('TR_%s_normal',   s.name); end
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
            e = base(c.px_per_deg, s.framerate, s.v_dps);
            e.experiment='TimeReversal_v25'; e.stim_type=s.name;
            e.reversed=rev; e.period_deg=s.period_deg;
            e.period_px=s.period_deg*c.px_per_deg;
            e.ratio_fly=s.ratio_fly;
            e.cycles_per_frame=(s.v_dps/s.period_deg)/s.framerate;
            e.gt_direction_deg = 0;
            dnn_v25_save(edir, dnn_v25_broadcast(mI, N_rows), e);
            count = count + 1;
        end
    end
    fprintf('  +TimeReversal: %d total (%d skipped)\n\n', count, skipped);
end
fprintf('Done.\n');

function s = polarity_label(iv)
% Human-readable polarity label, stored alongside the frames.
    if iv == 1
        s = 'barra oscura sobre fondo claro';
    else
        s = 'barra clara sobre fondo oscuro';
    end
end
