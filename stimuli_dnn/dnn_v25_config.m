function cfg = dnn_v25_config(experiment, N_cols)
% DNN_V25_CONFIG  Sampling parameters for the DNN exports (v25).
%
%   cfg = dnn_v25_config(experiment, N_cols)
%
%   Single source of truth for px_per_deg, framerate and clip length of
%   every stimulus family exported to the optical-flow networks.
%
%   DESIGN RULE
%   -----------
%   The stimulus itself (speed in deg/s, spatial frequency in cpd, period,
%   temporal envelope) is NEVER changed: those are the physical quantities
%   stINRF3 was run with. Only the SAMPLING is changed:
%
%       v_px_per_frame = speed_dps * px_per_deg / framerate
%
%   px_per_deg and framerate are display properties, so altering them
%   resamples the same physical stimulus rather than replacing it. They are
%   chosen so that
%
%       (a) 10 <= v_px_per_frame <= 40      the interframe displacement
%                                           range these networks were
%                                           trained on (MPI-Sintel: mean
%                                           13.5, p75 12.8, p90 39.0)
%       (b) cycles_per_frame <= 0.25        4x below Nyquist, so no
%                                           direction ambiguity survives
%                                           the networks' internal x8
%                                           downsampling
%       (c) n_frames >= 15                  Dual consumes 15-frame windows
%                                           (ME-Attention 11); below that
%                                           the multiframe models cannot
%                                           run at all
%
%   n_frames is set to 36 throughout (Dual gets 22 windows, ME-Attention
%   26), except DH, which keeps stINRF3's own 500 ms so that the Gaussian
%   temporal envelope occupies the same fraction of the clip in both.
%
% Inputs
%   experiment (char) 'gratings' | 'plaidI' | 'plaidII' | 'DH' |
%                     'Contrast' | 'Masking' | 'Missing' | 'Opponency' |
%                     'ReversePhi' | 'SecondOrder' | 'TimeReversal'
%   N_cols (int) canvas width in px, used to derive dov
%
% Outputs
%   cfg (struct) .px_per_deg .framerate .duration_ms .n_frames .dov
%                .notes, plus per-experiment grids

    cfg = struct();

    switch lower(experiment)

    % ---------------------------------------------------------------
    case 'gratings'
        cfg.px_per_deg  = 128/1.5;      % 85.333, unchanged
        cfg.framerate   = 40;           % was 200
        cfg.n_frames    = 36;           % 900 ms
        cfg.f_spatial   = 0.5;          % cpd, unchanged
        cfg.speed_Hz    = 8;            % unchanged -> 16 deg/s
        cfg.directions_deg = 0:15:90;
        cfg.notes = '34.1 px/frame, 0.20 cy/frame';

    % ---------------------------------------------------------------
    case 'texture'
        % Same physical stimulus as 'gratings' -- 0.5 cpd, 8 Hz, 85.333
        % px/deg, so still 16 deg/s -- resampled to 16 px/frame instead of
        % 34.13. Used by the four texture families (gratingTex, gratingMix,
        % gratingEdge, gratingSlit, barTex).
        %
        % WHY 16 AND NOT 34.13
        % --------------------
        % At 34.13 px/frame no texture serves all five networks, and this
        % was measured rather than assumed:
        %
        %   PWC-Net and LiteFlowNet2 are coarse-to-fine pyramids. They need
        %   structure finer than about 68 px to refine with; given only
        %   coarser texture they return ~100-200 px of flow, and do so even
        %   for two identical frames.
        %
        %   ME-Attention is a motion-energy model. A texture component of
        %   period lambda moving at v px/frame sits at v/lambda cycles per
        %   frame, and above 0.5 it is temporally aliased -- which for an
        %   energy model is not noise but energy pointing the wrong way. At
        %   v = 34.13 everything below 68 px aliases.
        %
        % So one needs components below 68 px and no components below 68 px
        % at the same time. Bands of 32, 68, 80 and 96 px were all tried:
        % 32 fixes the pyramids and breaks ME-Attention, the other three do
        % the reverse. There is no window.
        %
        % The aliasing threshold is 2*v, so lowering the speed opens one.
        % At 16 px/frame it falls to 32 px, and the usable band floor was
        % then measured directly, on pure translating texture at 45 deg:
        %
        %   floor  max temporal freq   result
        %    32 px      0.50 cy/fr     all five work, but the finest
        %                              component sits exactly ON Nyquist
        %    40 px      0.40 cy/fr     all five work, with margin
        %    48 px      0.33 cy/fr     PWC-Net -60.6, LiteFlowNet2 -78.1
        %    64 px      0.25 cy/fr     both still broken
        %
        % So the window is roughly 32 to 47 px and 40 sits in the middle.
        % At 40: RAFT 45.1, PWC-Net 45.0, LiteFlowNet2 45.0,
        % ME-Attention 42.3, Dual 43.6, against 45.0 true.
        %
        % 16 px/frame is still inside the 10-40 band, and cycles_per_frame
        % drops from 0.20 to 0.094, further from Nyquist than before.
        cfg.px_per_deg  = 128/1.5;      % 85.333, unchanged
        cfg.framerate   = 85.3333;      % was 40
        cfg.n_frames    = 36;           % 422 ms
        cfg.f_spatial   = 0.5;          % cpd, unchanged
        cfg.speed_Hz    = 8;            % unchanged -> 16 deg/s
        cfg.directions_deg = 0:15:90;
        cfg.notes = '16.0 px/frame, 0.094 cy/frame, texture band 40-60 px';
        % The upper end was 480 px at first, on the assumption that the
        % texture had to supply coarse structure of its own to fix the
        % correspondence. With the grating present it does not: the grating
        % is that coarse structure, at 170.7 px. Measured at 45 deg with
        % grating 0.4 + texture 0.12, narrowing the texture to 40-44,
        % 40-60 or 40-120 px changes nothing against the full 40-480
        % (RAFT 45.1/45.1/45.0/44.9, PWC-Net 44.9/44.8/44.9/44.9,
        % LiteFlowNet2 44.9 throughout, Dual 43.6/44.4/44.3/44.3,
        % ME-Attention 40.3/39.2/41.4/41.9, true 45.0). The narrow band is
        % kept because it says what the stimulus actually needs: fine
        % detail, and nothing else.

    % ---------------------------------------------------------------
    case 'plaidi'
        cfg.px_per_deg  = 128/1.5;
        cfg.framerate   = 40;           % was 200
        cfg.n_frames    = 36;
        cfg.f_spatial   = 2.0;          % cpd, unchanged
        cfg.speed_Hz    = 8;
        cfg.angle_offset = 45;
        cfg.rigid_directions_deg = 0:15:90;
        cfg.notes = '12.1 px/frame pattern, 8.5 component, 0.20 cy/frame';

    % ---------------------------------------------------------------
    case 'plaidii'
        cfg.px_per_deg  = 128/1.5;
        cfg.framerate   = 94;           % was 200
        cfg.n_frames    = 36;
        cfg.f_spatial   = 2.0;
        cfg.ioc_targets_deg = -30:5:30;
        cfg.speed_ratio = 2;
        cfg.comp_sep_deg = 20;
        cfg.families    = 'both';
        cfg.notes = '13.9 px/frame pattern, 8.5 fast component, 0.20 cy/frame';

    % ---------------------------------------------------------------
    case 'dh'
        % The only family whose px_per_deg has to move: with ppd = 100 the
        % 3 cpd moving component tops out at 8.3 px/frame and cannot reach
        % the band. ppd = 180 puts both components in it (12 and 36) at the
        % same 0.20 cy/frame, and brings dov to 4.27 deg, essentially the
        % dov = 4.0 that stINRF3_config declares for DH.
        cfg.px_per_deg  = 180;          % was 100
        cfg.framerate   = 60;           % was 1000
        cfg.duration_ms = 500;          % stINRF3's own duration -> 30 frames
        cfg.n_frames    = 30;
        cfg.tf_hz       = 12;           % see make_DH_stimulus_1D_v25 (bug B2)
        cfg.sfm         = [1, 3, 1, 3, 3];
        cfg.sfs         = [1, 1, 3, 1, 1];
        cfg.cond_names  = {'S1','S3','CS','CSr','CM'};
        cfg.sigma_t_ms  = [10, 25, 37.5, 100];
        cfg.notes = ['1 cpd -> 36 px/frame, 3 cpd -> 12 px/frame. ' ...
                     'sigma_t = 10 ms spans 2 frames: fine for the 2-frame ' ...
                     'models, only 2 of 15 frames carry signal for Dual.'];

    % ---------------------------------------------------------------
    case {'contrast','masking','secondorder'}
        % ONE framerate for the whole experiment. Using a different one per
        % spatial frequency would make px/frame co-vary with sf, so a
        % difference in performance across the sweep could not be
        % attributed to either axis. The finest sf sets the ceiling: at
        % ppd = 100 a 2 cpd grating has a 50 px period, so a quarter cycle
        % is 12.5 px -- and the band starts at 10. That leaves a narrow but
        % usable window, and 82 fps lands in it at 11.4 px/frame for every
        % condition (cy/frame 0.057 / 0.114 / 0.229).
        cfg.px_per_deg  = 100;          % unchanged
        cfg.speed_dps   = 9.38;         % unchanged
        cfg.sf_list     = [0.5 1 2];    % 4 cpd removed: ceiling 6.25 px/frame
        cfg.framerate   = 82;
        cfg.n_frames    = 36;
        if strcmpi(experiment, 'secondorder')
            % Only two spatial frequencies here, so the ceiling is looser
            % and the whole experiment fits at 20 px/frame.
            cfg.sf_list    = [0.5 1];
            % stINRF3 runs second order at 8 deg/s; 40 fps keeps that
            % physical speed at 20 px/frame. (Was 47 fps with the 9.38 deg/s
            % inherited from LightEdges, which is not this experiment's speed.)
            cfg.framerate  = 40;
            cfg.speed_dps  = 8;
            cfg.carrier_grain_px = 16;  % 2x the networks' x8 stride
            cfg.sigma_t_ms = 10000;
        elseif strcmpi(experiment, 'contrast')
            % stINRF3 runs contrast at 8 deg/s; 80 fps gives 10 px/frame.
            cfg.speed_dps  = 8;
            cfg.framerate  = 80;
            cfg.contrast_list = logspace(-2, 0, 20);   % the 20 contrasts of run_Contrast_1D.m
            cfg.sigma_t_ms = 100;
        else
            % stINRF3 runs masking at 10 deg/s; at 82 fps that is 12.2 px/frame.
            cfg.speed_dps   = 10;
            cfg.mask_ratios = [0.02 0.05 0.2 0.25 0.5 1 1.5 2 2.5 3 4 6 12];   % same 13 as run_Masking_1D.m
            cfg.sigma_t_ms = 1000;
            cfg.tf_mask_hz = 10;
            cfg.contrast_signal = 0.8;
            cfg.contrast_mask   = 0.8;
            % The mask drifts at 10 Hz; 82 fps puts it at 0.122 cy/frame.
        end
        cfg.notes = 'single framerate per experiment; sf = 4 cpd dropped';

    % ---------------------------------------------------------------
    case 'missing'
        % Stepped motion, 40 ms per step. Exporting at exactly 25 fps puts
        % one frame per step, so every consecutive pair shows exactly one
        % jump; at any other rate the 2-frame models get pairs with zero
        % displacement. sf is capped because the harmonic series runs to
        % 39*sf and must stay below the pixel Nyquist (ppd/2 = 50 cpd).
        cfg.px_per_deg  = 100;          % unchanged
        cfg.framerate   = 25;           % = 1/0.040 s
        cfg.speed_dps   = 6.0;          % stINRF3's speed -> 24 px per 40 ms step
        cfg.sf_list     = [0.5 1.0];    % 2 and 4 cpd alias spatially
        cfg.n_frames    = 36;           % 1440 ms
        cfg.types       = [1, -1];
        cfg.sigma_t_ms  = 10000;
        cfg.notes = '20 px per 40 ms step; 39*sf < 50 cpd Nyquist';

    % ---------------------------------------------------------------
    case 'opponency'
        % ppd 100 -> 200 so the narrowest bar (0.125 deg) is 25 px and can
        % take a 12 px displacement; below that the bar jumps more than
        % half its own width and correspondence breaks. One ppd for all
        % four widths keeps them comparable. dov = 3.84 deg, close to
        % stINRF3's 4.0.
        cfg.px_per_deg  = 200;          % was 100
        cfg.framerate   = 400/3;        % 133.33 fps (was 156): 12 px/frame at 8 deg/s
        cfg.speed_dps   = 8;            % stINRF3's speed (was the inherited 9.38)
        cfg.bar_widths  = [0.125 0.25 0.5 1.0];
        cfg.directions  = [1, -1];
        % run_Opponency_1D.m recorre las dos polaridades de contraste:
        % inver = 1 -> barra oscura sobre fondo claro
        % inver = 0 -> barra clara sobre fondo oscuro
        cfg.inver       = [1 0];
        cfg.n_frames    = 36;
        cfg.sigma_t_ms  = 100000;
        cfg.notes = '12 px/frame for all four widths (limited by 0.125 deg)';

    % ---------------------------------------------------------------
    case 'reversephi'
        % The effect lives in the frame-to-frame contrast inversion, so the
        % clip cannot be decimated: it is generated at the low rate
        % directly. White noise of 1 px grain has no correspondence at all
        % across a 20 px jump, hence grain >= 2x displacement.
        cfg.px_per_deg  = 25;           % unchanged
        cfg.framerate   = 11.7;
        cfg.speed_dps   = 9.38;         % unchanged
        cfg.grain_px    = 50;           % was 1 px
        cfg.contrast_list = [0.25 0.5 0.75 1.0];
        cfg.types       = [1, -1];
        cfg.n_frames    = 36;           % 3077 ms
        cfg.sigma_t_ms  = 10000;
        cfg.notes = '20 px/frame, grain 50 px = 2x displacement';

    % ---------------------------------------------------------------
    case 'timereversal'
        % Both stimuli at stINRF3's physical speed and period, 10 px/frame.
        % LightEdges: 9.375 deg/s (stINRF3 9.38) at 170.67 ppd and 160 fps.
        % UnevenSquareWave: 1.17 deg/s at 85.33 ppd needs 9.984 fps; 15
        % frames (1.5 s) so that Dual gets its 15-frame input window.
        cfg.stim = struct( ...
            'name',       {'LightEdges',        'UnevenSquareWave'}, ...
            'px_per_deg', {512/3,               128/1.5}, ...
            'v_dps',      {9.375,               1.17}, ...
            'period_deg', {1.17,                2.34}, ...
            'framerate',  {160,                 1.17*(128/1.5)/10}, ...
            'n_frames',   {36,                  15}, ...
            'ratio_fly',  {5.72,                12.74});
        cfg.px_per_deg = 128/1.5;       % default; each stimulus carries its own
        cfg.n_frames   = 36;
        cfg.notes = '10 px/frame both, physical speed as in stINRF3';

    otherwise
        error('dnn_v25_config: unknown experiment ''%s''', experiment);
    end

    % ---- derived --------------------------------------------------
    if nargin >= 2 && ~isempty(N_cols)
        cfg.N_cols = N_cols;
        cfg.dov    = N_cols / cfg.px_per_deg;
    end
    if ~isfield(cfg, 'duration_ms')
        if isfield(cfg, 'framerate')
            cfg.duration_ms = 1000 * cfg.n_frames / cfg.framerate;
        elseif isfield(cfg, 'framerates')
            % one framerate per spatial frequency -> one duration each
            cfg.durations_ms = 1000 * cfg.n_frames ./ cfg.framerates;
            cfg.duration_ms  = max(cfg.durations_ms);
        elseif isfield(cfg, 'stim')
            frs = [cfg.stim.framerate];
            for k = 1:numel(cfg.stim)
                cfg.stim(k).duration_ms = 1000 * cfg.n_frames / frs(k);
            end
            cfg.duration_ms = max(1000 * cfg.n_frames ./ frs);
        end
    end
end
