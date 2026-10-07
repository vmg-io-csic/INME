%% export_dnn_v25_gratingtex.m
%
% Drifting grating with a rigidly co-moving texture added on top, for the
% five optical-flow networks.
%
% WHY
% ---
% The companion experiment (export_dnn_v25_bartex.m) changes two things at
% once relative to the bare grating: it adds trackable 2D structure AND it
% adds a figure/ground boundary. If the networks succeed there, that alone
% does not say which of the two they needed.
%
% This one changes only the first. The grating stays exactly as it is, and
% a band-limited texture is added that translates with precisely the same
% velocity vector, so the composite still has a single global motion and
% still has no object boundary anywhere -- but it now carries features that
% can be tracked, which the bare sinusoid does not. The aperture problem is
% removed and nothing else is.
%
% The texture contrast is swept, and c_t = 0 reproduces the bare grating
% bit for bit, so the ladder carries its own control: whatever the networks
% do at c_t = 0 is what they already do on the grating stimuli, and any
% change along the ladder is attributable to the texture alone.
%
% SAMPLING
% --------
% dnn_v25_config('texture'): 85.333 px/deg, 85.33 fps, 36 frames,
% 16.00 px/frame, 0.5 cpd, directions 0:15:90. Unchanged.
%
% DESIGN CHOICES
% --------------
%   grating contrast 0.4
%       Held fixed across the ladder. Room has to be left for the texture,
%       and the grating's contrast must not co-vary with it or the two
%       axes could not be told apart.
%   texture contrast 0, 0.04, 0.08, 0.12 (RMS)
%       Over 36 frames the grating drifts through every phase relative to
%       the static texture pattern, so the two peaks do eventually
%       coincide and the worst case is the sum of both: the band-limited
%       field peaks at 3.82 sigma, so the constraint is
%       c_grating + peak*c_texture <= 1. The peak of a unit-RMS
%       band-limited field is not a universal constant: it depends on the
%       canvas size and the seed, and over the three canvases it runs from
%       3.92 to 4.65 sigma. Taking the worst case, 4.65, the ladder has to
%       stop at 0.4 + 4.65*0.12 = 0.958. Every stimulus is checked against
%       info.clipped as it is written, so a bad constant fails loudly
%       rather than silently flattening the peaks. Every
%       stimulus is checked against info.clipped as it is written.
%   texture band 40-60 px per cycle
%       Bracketing the grating's own 170.7 px period at the top. The bottom
%       was 96 px at first, so that every component stayed above twice the
%       16.00 px displacement, but that starved PWC-Net and LiteFlowNet2 of
%       the near-pixel detail their coarse-to-fine pyramids refine with:
%       they returned ~100 px of flow even for two identical frames.
%       Extending to 32 px fixes both and leaves RAFT unchanged. A 32-96 px
%       component is ambiguous alone at this displacement; the composite is
%       not, because the coarse component fixes the correspondence.
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

meanL       = 100;
c_grating   = 0.4;
c_tex_list  = [0 0.04 0.08 0.12];
period_lo   = 40;
period_hi   = 60;
seed_tex    = 2;
phase0      = 0;

for g = 1:numel(groups)
    tag = groups(g).tag;  N_cols = groups(g).N_cols;  N_rows = groups(g).N_rows;
    out_root = fullfile('stimuli_by_model', tag, 'v25');
    c = dnn_v25_config('texture', N_cols);
    f_cpx = c.f_spatial / c.px_per_deg;
    speed = c.speed_Hz / c.framerate;
    v_px  = speed / f_cpx;
    fprintf('=== %s (%dx%d)  v = %.3f px/frame ===\n', tag, N_cols, N_rows, v_px);

    count = 0; skipped = 0;
    for ct = c_tex_list
        for d = c.directions_deg
            name = sprintf('gratingTex_dir%03d_ct%03d', d, round(ct*1000));
            edir = fullfile(out_root, name);
            if exist(fullfile(edir,'frames.mat'),'file')
                skipped = skipped + 1;  count = count + 1;  continue;
            end
            [TCB, info] = make_textured_grating(N_cols, c.n_frames, f_cpx, ...
                    speed, d, phase0, c_grating, ct, meanL, ...
                    period_lo, period_hi, seed_tex);
            if info.clipped
                error('gratingTex %s clipped: lower c_texture', name);
            end
            extra = struct('experiment','gratingTex_v25', ...
                'direction_deg', d, 'gt_direction_deg', d, ...
                'grating_contrast', c_grating, ...
                'texture_rms_contrast', ct, ...
                'texture_period_lo_px', period_lo, ...
                'texture_period_hi_px', period_hi, ...
                'seed_texture', seed_tex, ...
                'f_spatial_cpd', c.f_spatial, 'speed_Hz', c.speed_Hz, ...
                'px_per_deg', c.px_per_deg, 'framerate', c.framerate, ...
                'dov', c.dov, 'N_cols', N_cols, 'N_rows', N_rows, ...
                'meanL', meanL, ...
                'speed_px_per_frame', v_px, ...
                'cycles_per_frame', c.speed_Hz / c.framerate, ...
                'max_excursion', info.max_excursion);
            dnn_v25_save(edir, dnn_v25_crop(TCB, N_rows), extra);
            count = count + 1;
        end
    end
    fprintf('  gratingTex: %d (%d ya existian)\n', count, skipped);
end
fprintf('Done.\n');
