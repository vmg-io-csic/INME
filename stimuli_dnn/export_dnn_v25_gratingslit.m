%% export_dnn_v25_gratingslit.m
%
% Drifting grating seen through a slit that runs ALONG the motion, so the
% bars are cut across and their end points are visible.
%
% This is the condition export_dnn_v25_gratingedge.m was meant to be. That
% file puts the aperture across the motion, which makes its edges parallel
% to the bars: nothing gets cut, no end point exists, and the stimulus adds
% a border without adding any two-dimensional information. It is kept, but
% as the control it actually is -- border without terminators -- and this
% file supplies the informative case.
%
%   The pair is the point:
%     gratingEdge   aperture across the motion, edges parallel to the bars
%                   -> a border, no terminators
%     gratingSlit   aperture along the motion, edges across the bars
%                   -> a border WITH terminators moving at the true velocity
%
%   A model that recovers the direction from gratingSlit but not from
%   gratingEdge is using the terminators, not the mere presence of an edge.
%
% Sampling from dnn_v25_config('texture') and
% grating contrast 0.4, the same as gratingTex_ct000 and gratingEdge, so
% the three differ only in what the aperture does.
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
c_grating  = 0.4;
slit_px    = [192 384];
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
    for w = slit_px
        for d = c.directions_deg
            name = sprintf('gratingSlit_dir%03d_w%03d', d, w);
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file')
                skipped = skipped + 1;  count = count + 1;  continue;
            end
            [TCB, info] = make_slit_grating(N_cols, c.n_frames, f_cpx, ...
                    speed, d, phase0, c_grating, meanL, w);
            if info.clipped, error('gratingSlit %s clipped', name); end
            extra = struct('experiment','gratingSlit_v25', ...
                'direction_deg', d, 'gt_direction_deg', d, ...
                'slit_width_px', w, 'slit_width_deg', w/c.px_per_deg, ...
                'grating_contrast', c_grating, ...
                'duty_cycle', info.duty_cycle, ...
                'f_spatial_cpd', c.f_spatial, 'speed_Hz', c.speed_Hz, ...
                'px_per_deg', c.px_per_deg, 'framerate', c.framerate, ...
                'dov', c.dov, 'N_cols', N_cols, 'N_rows', N_rows, ...
                'meanL', meanL, 'speed_px_per_frame', v_px, ...
                'cycles_per_frame', c.speed_Hz / c.framerate);
            dnn_v25_save(edir, dnn_v25_crop(TCB, N_rows), extra);
            count = count + 1;
        end
    end
    fprintf('  gratingSlit: %d (%d ya existian)\n', count, skipped);
end
fprintf('Done.\n');
