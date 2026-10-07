%% export_dnn_v25_gratings.m
%
% Exports gratings, Type-I plaids and Type-II plaids for the five
% optical-flow networks, resampled into the 10-40 px/frame band.
%
% NOTHING about the stimulus changes: spatial frequency stays at 0.5 cpd
% (gratings) and 2.0 cpd (plaids), the temporal frequency stays at 8 Hz,
% px_per_deg stays at 85.333. Only the framerate is lowered, which is
% resampling, not redesign.
%
%   framerate = 5 * tf_max_Hz   ->  exactly 0.20 cycles/frame
%
%   gratings  8 Hz  -> 40 fps -> 34.1 px/frame
%   plaidI    8 Hz  -> 40 fps -> 12.1 px/frame (pattern), 8.5 (component)
%   plaidII  18.75Hz-> 94 fps -> 13.9 px/frame (pattern), 8.5 (fast comp)
%
% Clips are 36 frames so Dual (15-frame windows) gets 22 windows and
% ME-Attention (11) gets 26.
%
% Resume-safe: skips any stimulus whose frames.mat already exists.
% Set SWEEP_ONLY_TAG=<tag> to process one canvas group only.

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

meanL = 100;  contrast = 1;  phase0 = 0;

for g = 1:numel(groups)
    tag = groups(g).tag;  N_cols = groups(g).N_cols;  N_rows = groups(g).N_rows;
    out_root = fullfile('stimuli_by_model', tag, 'v25');
    fprintf('=== %s (%dx%d) ===\n', tag, N_cols, N_rows);
    count = 0; skipped = 0;

    %% ---------------- GRATINGS ----------------
    c = dnn_v25_config('gratings', N_cols);
    f_cpx = c.f_spatial / c.px_per_deg;
    speed = c.speed_Hz / c.framerate;            % cycles/frame
    for d = c.directions_deg
        name = sprintf('gratings_dir%03d', d);
        edir = fullfile(out_root, name);
        if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
        TCB = make_grating(N_cols, c.n_frames, f_cpx, speed, d, phase0, contrast, meanL);
        extra = struct('experiment','gratings_v25','direction_deg',d, ...
            'f_spatial_cpd',c.f_spatial,'speed_Hz',c.speed_Hz, ...
            'px_per_deg',c.px_per_deg,'framerate',c.framerate,'dov',c.dov, ...
            'N_cols',N_cols,'N_rows',N_rows,'meanL',meanL, ...
            'speed_px_per_frame', c.px_per_deg*(c.speed_Hz/c.f_spatial)/c.framerate, ...
            'cycles_per_frame', c.speed_Hz/c.framerate, ...
            'gt_direction_deg', d);
        dnn_v25_save(edir, dnn_v25_crop(TCB, N_rows), extra);
        count = count + 1;
    end
    fprintf('  gratings: %d (%d skipped)\n', count, skipped);

    %% ---------------- PLAID I ----------------
    c = dnn_v25_config('plaidI', N_cols);
    f_cpx = c.f_spatial / c.px_per_deg;
    speed = c.speed_Hz / c.framerate;
    v_comp_px = c.px_per_deg*(c.speed_Hz/c.f_spatial)/c.framerate;
    for d = c.rigid_directions_deg
        th1 = d - c.angle_offset;  th2 = d + c.angle_offset;
        name = sprintf('plaidI_dir%03d', d);
        edir = fullfile(out_root, name);
        if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
        TCB = make_plaidI(N_cols, c.n_frames, f_cpx, speed, th1, th2, ...
                          phase0, phase0, contrast, meanL);
        extra = struct('experiment','plaidI_v25','rigid_dir_deg',d, ...
            'theta1_deg',th1,'theta2_deg',th2, ...
            'f_spatial_cpd',c.f_spatial,'speed_Hz',c.speed_Hz, ...
            'px_per_deg',c.px_per_deg,'framerate',c.framerate,'dov',c.dov, ...
            'N_cols',N_cols,'N_rows',N_rows,'meanL',meanL, ...
            'speed_px_per_frame', v_comp_px/cosd(c.angle_offset), ...
            'component_px_per_frame', v_comp_px, ...
            'cycles_per_frame', c.speed_Hz/c.framerate, ...
            'gt_direction_deg', d);
        dnn_v25_save(edir, dnn_v25_crop(TCB, N_rows), extra);
        count = count + 1;
    end
    fprintf('  plaidI: %d cumulative\n', count);

    %% ---------------- PLAID II ----------------
    c = dnn_v25_config('plaidII', N_cols);
    f_cpx = c.f_spatial / c.px_per_deg;
    % tf_scale converts the table's relative speeds to cycles/frame; the
    % 0.21333 reference is the original grating speed (see run_plaidII_2D).
    tf_scale = c.f_spatial / 0.21333 / c.framerate;
    [configs, geo] = make_plaidII_configs(c.ioc_targets_deg, c.speed_ratio, ...
                                          c.comp_sep_deg, c.families);
    for k = 1:size(configs,1)
        s1 = configs(k,1)*tf_scale;  s2 = configs(k,2)*tf_scale;
        th1 = configs(k,3);  th2 = configs(k,4);  ioc = configs(k,5);
        name = sprintf('plaidII_ioc%+04d_c%02d', round(ioc), k);
        edir = fullfile(out_root, name);
        if exist(fullfile(edir,'frames.mat'),'file'), skipped=skipped+1; count=count+1; continue; end
        TCB = make_plaidII(N_cols, c.n_frames, f_cpx, s1, s2, th1, th2, ...
                           phase0, phase0, contrast, meanL);
        % pattern speed from the IOC geometry: |V| = v_i / cos(delta_i)
        d1 = th1 - ioc;
        v1_px = s1 * c.px_per_deg / c.f_spatial;      % px/frame, component 1
        extra = struct('experiment','plaidII_v25','ioc_dir_deg',ioc, ...
            'theta1_deg',th1,'theta2_deg',th2,'speed1',configs(k,1), ...
            'speed2',configs(k,2),'speed_ratio',c.speed_ratio, ...
            'comp_sep_deg',c.comp_sep_deg, ...
            'f_spatial_cpd',c.f_spatial, ...
            'px_per_deg',c.px_per_deg,'framerate',c.framerate,'dov',c.dov, ...
            'N_cols',N_cols,'N_rows',N_rows,'meanL',meanL, ...
            'speed_px_per_frame', v1_px/cosd(d1), ...
            'component_px_per_frame', max(s1,s2)*c.px_per_deg/c.f_spatial, ...
            'cycles_per_frame', max(s1,s2), ...
            'gt_direction_deg', ioc);
        dnn_v25_save(edir, dnn_v25_crop(TCB, N_rows), extra);
        count = count + 1;
    end
    fprintf('  plaidII: %d total (%d skipped)\n\n', count, skipped);
end
fprintf('Done.\n');
