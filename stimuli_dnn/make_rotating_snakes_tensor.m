function [TCB, I_lum, info] = make_rotating_snakes_tensor(N, period_px, ...
                                  n_frames, onset_f, offset_f, meanL, opts)
% MAKE_ROTATING_SNAKES_TENSOR  Synthetic Rotating-Snakes stimulus tensor.
%
%   [TCB, I_lum, info] = make_rotating_snakes_tensor(N, period_px, ...
%                            n_frames, onset_f, offset_f, meanL, opts)
%
%   Builds the Akiyoshi Kitaoka "Rotating Snakes" pattern PARAMETRICALLY
%   (single concentric wheel centred on the canvas) instead of downloading
%   the copyrighted original image. The asymmetric snake unit is a 4-step
%   simplified luminance profile repeated angularly; concentric rings
%   alternate the gradient direction so adjacent rings appear to rotate in
%   opposite senses - the canonical Rotating-Snakes signature.
%
%   NO STIMULUS-SIDE BLUR. All spatial smoothing is delegated to the
%   model's A-stage (sigmaMu in stINRF3_compute). To compensate for the
%   removed stimulus blur, set model.sigmaMu = 0.12 deg in the
%   rotating_snakes case of stINRF3_config.m.
%
%   RING LAYOUT - exact integer fit. Each ring must fit a whole number of
%   snake units around its circumference. A naive round() of the unit
%   count leaves a residual detuning (the pattern is slightly compressed
%   or stretched vs period_px) whose sign alternates across rings; stINRF3,
%   being frequency-tuned, reads that as a spurious alternating signal.
%   This version chooses INTEGER unit counts first and derives each ring's
%   exact radius from them (r = n*period_px/(2*pi)), so the detuning is
%   exactly zero in every ring.
%
%   Inputs
%   ------
%   N          [int]     Canvas size (square). Use 300 to match stINRF3
%                        calling convention.
%   period_px  [double]  Tangential period of one snake unit [px]
%                        (= dpi/sf0_cpd, the period stINRF3 is tuned to).
%   n_frames   [int]     Total number of frames.
%   onset_f    [int]     First visible frame.
%   offset_f   [int]     Last visible frame.
%   meanL      [double]  Background / mean luminance (e.g. 50).
%   opts       [struct]  Optional, all fields optional:
%       .n_rings        number of concentric snake rings        (default 3)
%       .r_inner        inner radius of ring stack [px]    (default 0.15*N)
%       .r_outer        outer radius of ring stack [px]    (default 0.46*N)
%       .levels         4-vector, the simplified luminance steps in [0,100]
%                                                  (default [0 25 100 50])
%       .gap_frac       blank-grey gap between rings, fraction of ring width
%                                                           (default 0.10)
%       .uniform_sense  if true, all rings share gradient sense +1
%                       (diagnostic control). Default: false.
%
%   Outputs
%   -------
%   TCB    [single, N x N x n_frames]  Stimulus tensor for stINRF3_compute.
%   I_lum  [N x N double]              The static pattern in [0,100].
%   info   struct                      .period_px .n_rings .levels
%                                       .ring_edges .center .gradient_sense
%                                       .n_units .round_error
%
%   Notes
%   -----
%   - r_inner / r_outer are TARGETS; the actual ring radii snap to exact
%     integer-unit values, so they may shift slightly.
%   - info.n_units lists the exact unit count per ring; info.round_error
%     is ~0 by construction.
%   - Pass the SAME N to stINRF3_compute as used here.

    if nargin < 7 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'n_rings'),       opts.n_rings        = 3;            end
    if ~isfield(opts,'r_inner'),       opts.r_inner        = 0.15 * N;     end
    if ~isfield(opts,'r_outer'),       opts.r_outer        = 0.46 * N;     end
    if ~isfield(opts,'levels'),        opts.levels         = [0 25 100 50]; end
    if ~isfield(opts,'gap_frac'),      opts.gap_frac       = 0.10;         end
    if ~isfield(opts,'uniform_sense'), opts.uniform_sense  = false;        end
    if ~isfield(opts,'angular_offset'), opts.angular_offset = 0;           end
    if ~isfield(opts,'n_units_fixed'), opts.n_units_fixed  = [];           end
    % uniform_sense:
    %   false / 0       -> alternating  [+1 -1 +1 ...] (canonical Kitaoka)
    %   true  / +1      -> all +1       (uniform CCW polarity)
    %         / -1      -> all -1       (uniform CW polarity)
    %
    % angular_offset: progressive azimuthal shift between consecutive rings,
    %   measured as a fraction of one snake unit. 0 = rings synchronised
    %   (default). 0.5 = each ring shifted half a unit relative to the
    %   previous one (reproduces the curved "spiral" appearance of
    %   Kitaoka 2003 Fig. 4).
    %
    % n_units_fixed: if empty (default), n_units per ring is derived from
    %   period_px so that the tangential period at each ring's radius
    %   matches period_px exactly. If set to a positive integer, all rings
    %   use that fixed number of angular units regardless of radius. With
    %   fixed n_units the effective period_px varies with the radius
    %   (shorter in inner rings, longer in outer rings), so most rings
    %   will be out of the model's spatial-frequency tuning. Use only
    %   when reproducing Kitaoka-style geometries with many concentric
    %   rings sharing the same angular layout (Fig. 4).

    levels = opts.levels(:).';
    if numel(levels) ~= 4
        error('make_rotating_snakes_tensor:levels', ...
              'opts.levels must be a 4-vector (got %d elements).', numel(levels));
    end

    % ---- canvas geometry (symmetric centre) ------------------------------
    cx = (N + 1) / 2;
    cy = (N + 1) / 2;
    [Xg, Yg] = meshgrid(1:N, 1:N);
    dx = Xg - cx;
    dy = Yg - cy;
    r  = sqrt(dx.^2 + dy.^2);
    th = atan2(dy, dx);                 % [-pi, pi], azimuth

    % ---- concentric ring layout: two modes -------------------------------
    %  (a) n_units_fixed empty  -> derive n_units from period_px per ring
    %                              (anchored to model's spatial tuning)
    %  (b) n_units_fixed set    -> all rings use the same n_units;
    %                              radii distributed evenly across the band
    n_rings   = opts.n_rings;
    edges0    = linspace(opts.r_inner, opts.r_outer, n_rings + 1);
    ring_w    = edges0(2) - edges0(1);
    gap       = opts.gap_frac * ring_w;
    half_band = 0.5 * (ring_w - gap);

    if isempty(opts.n_units_fixed)
        % Mode (a): n_units derived from period_px so each ring is at the
        % model's tuned spatial frequency. Radii snap to exact integer-unit
        % values.
        r_of_n = @(n) n * period_px / (2*pi);
        n_units_per_ring = zeros(1, n_rings);
        r_mid_per_ring   = zeros(1, n_rings);
        n_prev = 0;
        for k = 1:n_rings
            if k == 1
                r_need = opts.r_inner + half_band;
            else
                r_need = r_mid_per_ring(k-1) + 2*half_band + gap;
            end
            n_try = max(n_prev + 1, ceil(2*pi*r_need / period_px));
            n_units_per_ring(k) = n_try;
            r_mid_per_ring(k)   = r_of_n(n_try);
            n_prev = n_try;
        end
        round_error = n_units_per_ring - (2*pi*r_mid_per_ring / period_px);
    else
        % Mode (b): all rings share the same n_units (Kitaoka Fig. 4 style).
        % Radii distributed evenly between r_inner and r_outer. Effective
        % period varies with radius; most rings will be out of tuning.
        n_fixed = max(2, round(opts.n_units_fixed));
        n_units_per_ring = n_fixed * ones(1, n_rings);
        r_mid_per_ring = zeros(1, n_rings);
        for k = 1:n_rings
            r_mid_per_ring(k) = opts.r_inner + half_band + (k-1) * ring_w;
        end
        round_error = zeros(1, n_rings);   % no detuning to report in this mode
    end

    edges = zeros(1, n_rings + 1);
    edges(1) = r_mid_per_ring(1) - half_band - gap/2;
    for k = 1:n_rings
        edges(k+1) = r_mid_per_ring(k) + half_band + gap/2;
    end

    I_lum     = meanL * ones(N, N);
    grad_sense = zeros(1, n_rings);

    for k = 1:n_rings
        n_units = n_units_per_ring(k);
        r_mid   = r_mid_per_ring(k);
        r_lo    = r_mid - half_band;
        r_hi    = r_mid + half_band;
        ring_mask = (r >= r_lo) & (r < r_hi);

        % phase within one unit, [0,1) - n_units exact, no round error.
        % Progressive azimuthal offset applied per ring: each ring is
        % shifted by (k-1) * angular_offset units relative to ring 1.
        u = mod(th / (2*pi) * n_units - (k-1) * opts.angular_offset, 1);

        % gradient sense:
        %   alternating  (canonical Kitaoka wheel, opts.uniform_sense = false / 0)
        %   uniform +1   (CCW polarity, opts.uniform_sense = true / +1)
        %   uniform -1   (CW  polarity, opts.uniform_sense = -1)
        if islogical(opts.uniform_sense) || opts.uniform_sense == 0
            if opts.uniform_sense
                s = 1;
            else
                s = (-1)^(k+1);
            end
        else
            s = sign(opts.uniform_sense);   % +1 or -1
        end
        % Note on sign convention: the model's Gabor is built with y
        % decreasing along the row index (see stINRF3_compute), so the
        % sign of the model's response is inverted relative to the
        % stimulus's +Yg convention. By calibration against the
        % real-motion control (uniform tangential with v>0 produces
        % score<0 in the readout), the polarity stored in grad_sense
        % must be the OPPOSITE of the visually-painted polarity s.
        grad_sense(k) = -s;
        if s < 0
            u = 1 - u;                  % flip the asymmetric profile
        end

        % ---- simplified stepped luminance profile (4 steps) -------------
        % Kitaoka & Ashida 2003, Fig. 2: the spatial order in the direction
        % of perceived motion is  black -> dark-gray -> white -> light-gray.
        % With levels = [black, dark-gray, white, light-gray] = [0, 25, 100, 50],
        % we paint them in this order as u increases:
        ring_val = zeros(size(u));
        ring_val(u <  0.25)               = levels(1);   % black
        ring_val(u >= 0.25 & u < 0.50)    = levels(2);   % dark-gray
        ring_val(u >= 0.50 & u < 0.75)    = levels(3);   % white
        ring_val(u >= 0.75)               = levels(4);   % light-gray

        I_lum(ring_mask) = ring_val(ring_mask);
    end

    % NO blur: all spatial smoothing happens at the A-stage of
    % stINRF3_compute (sigmaMu = 0.12 deg for the rotating_snakes case).

    % ---- onset-paradigm tensor (grey -> pattern -> grey) -----------------
    TCB   = single(meanL * ones(N, N, n_frames));
    patch = single(I_lum);
    for t = onset_f : min(offset_f, n_frames)
        TCB(:,:,t) = patch;
    end

    info = struct('period_px', period_px, 'n_rings', n_rings, ...
                  'levels', levels, 'ring_edges', edges, ...
                  'center', [cx cy], 'gradient_sense', grad_sense, ...
                  'n_units', n_units_per_ring, ...
                  'angular_offset', opts.angular_offset, ...
                  'round_error', round_error);
end
