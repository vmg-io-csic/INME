%% export_dnn_v25_bartex.m
%
% Textured bar over a differently-textured background, for the five
% optical-flow networks.
%
% WHY
% ---
% On the bare sinusoidal gratings the networks return the correct speed
% normal to the bars and an unconstrained slide along them -- the aperture
% problem. stINRF3 recovers the direction from the same stimuli, so the
% failure is the networks', but the aperture problem explains its form and
% leaves open what the networks would need in order to succeed. This
% stimulus supplies the thing their training data always has and the
% grating never does: two regions of different texture separated by a
% moving object boundary.
%
% SAMPLING
% --------
% From dnn_v25_config('texture'): 85.333 px/deg, 85.33 fps, 36 frames,
% 16.00 px/frame, same seven directions. That
% is deliberate -- the only difference between this run and the grating run
% is the stimulus itself, so the comparison isolates it.
%
% DESIGN CHOICES
% --------------
%   texture band 32-480 px per cycle
%       The upper end brackets the grating's own 170.7 px period. The lower
%       end was 96 px at first, chosen so that every component stayed above
%       twice the 16.00 px displacement and was therefore unambiguous on
%       its own. That turned out to break PWC-Net and LiteFlowNet2: they
%       are coarse-to-fine pyramid networks that resolve large motion at
%       the low levels and refine with near-pixel detail, and with nothing
%       below 96 px they have nothing to refine with -- they returned about
%       100 px of flow even for two identical frames. RAFT, which uses
%       all-pairs correlation instead of a warping pyramid, was unaffected.
%
%       Extending the band down to 32 px fixes both (they recover
%       34.16-34.22 px/frame against 34.133 true). A 32-96 px component is
%       ambiguous in isolation at this displacement, but the stimulus is
%       not: the coarse component fixes the correspondence and the fine one
%       only refines it, which is exactly how those pyramids work. 32 px
%       also stays above the ~16 px that survives the x8 downsampling.
%   contrast 0.20
%       The peak of a unit-RMS band-limited field depends on canvas size
%       and seed; across the three canvases and the two seeds it runs from
%       3.92 to 4.65 sigma. Against the worst case, 0.20 puts the extremes
%       at 0.93 of full range. Every stimulus is checked against
%       info.clipped as it is written.
%   bar widths 192 and 384 px (2.25 and 4.50 deg)
%       Both at least twice the coarsest texture period, so each bar
%       actually contains texture rather than a fragment of one blob.
%   static background
%       The figure carries the only motion, as in the pilot that motivated
%       this experiment.
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

meanL        = 100;
contrast     = 0.20;
period_lo    = 40;
period_hi    = 480;
seed_bar     = 2;
seed_bg      = 1;
bg_v_px      = 0;            % static background
bar_widths_px = [192 384];

for g = 1:numel(groups)
    tag = groups(g).tag;  N_cols = groups(g).N_cols;  N_rows = groups(g).N_rows;
    out_root = fullfile('stimuli_by_model', tag, 'v25');
    c = dnn_v25_config('texture', N_cols);
    v_px = c.px_per_deg * (c.speed_Hz / c.f_spatial) / c.framerate;
    fprintf('=== %s (%dx%d)  v = %.3f px/frame ===\n', tag, N_cols, N_rows, v_px);

    count = 0; skipped = 0;
    for w = bar_widths_px
        for d = c.directions_deg
            name = sprintf('barTex_dir%03d_w%03d', d, w);
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file')
                skipped = skipped + 1;  count = count + 1;  continue;
            end
            [TCB, info] = make_textured_bar(N_cols, c.n_frames, v_px, d, w, ...
                    contrast, meanL, period_lo, period_hi, ...
                    seed_bar, seed_bg, bg_v_px);
            if info.clipped
                error('barTex %s clipped: lower the contrast', name);
            end
            extra = struct('experiment','barTex_v25', ...
                'direction_deg', d, 'gt_direction_deg', d, ...
                'bar_width_px', w, 'bar_width_deg', w/c.px_per_deg, ...
                'contrast', contrast, ...
                'texture_period_lo_px', period_lo, ...
                'texture_period_hi_px', period_hi, ...
                'seed_bar', seed_bar, 'seed_bg', seed_bg, ...
                'bg_v_px_per_frame', bg_v_px, ...
                'duty_cycle', info.duty_cycle, ...
                'px_per_deg', c.px_per_deg, 'framerate', c.framerate, ...
                'dov', c.dov, 'N_cols', N_cols, 'N_rows', N_rows, ...
                'meanL', meanL, ...
                'speed_px_per_frame', v_px, ...
                'wrap_period_px', N_cols);
            dnn_v25_save(edir, dnn_v25_crop(TCB, N_rows), extra);
            count = count + 1;
        end
    end
    fprintf('  barTex: %d (%d ya existian)\n', count, skipped);
end
fprintf('Done.\n');
