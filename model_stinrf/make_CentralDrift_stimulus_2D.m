function I_lum = make_CentralDrift_stimulus_2D(N, contrast, meanL, n_segments_override)
%MAKE_CENTRALDRIFT_STIMULUS_2D  Static central-drift illusion pattern.
%
%   Builds a ring of identical elements centred on the canvas. Each element
%   carries the luminance template returned by central_drift_template, whose
%   luminance runs as a continuous asymmetric gradient from dark to light;
%   repeating that gradient around the ring is what elicits illusory motion.
%   The pattern is entirely static: nothing moves between frames.
%
%   For every pixel the generator finds which angular segment it falls in,
%   expresses its position in the local tangential and radial frame of that
%   segment, and reads the template at that position by bilinear
%   interpolation, filling the background outside the elements.
%
%   The gradient runs in the direction of increasing polar angle, which reads
%   as clockwise. For a counter-clockwise condition, flip the returned image
%   at the call site rather than regenerating it with a reversed gradient, so
%   that the bias-cancelling readout used elsewhere still applies.
%
%   Inputs
%     N                    int. Canvas side in pixels; the canvas is square.
%     contrast             double. Scales the template's peak deviation from
%                          the background. contrast = 1 reproduces the
%                          template's own contrast.
%     meanL                double. Background luminance, which is also the
%                          luminance at the trailing edge of each element.
%     n_segments_override  int, optional. Place a different number of elements
%                          around the ring instead of the template's own
%                          count. Useful because the model's tangential
%                          readout cancels exactly for any arrangement that is
%                          invariant under a 90 degree rotation, so a multiple
%                          of four gives a degenerate, exactly zero response.
%
%   Outputs
%     I_lum  double, N x N. The static pattern, with the background at meanL.
%
%   Example
%     I = make_CentralDrift_stimulus_2D(512, 1, 100, 11);

    d = central_drift_template();

    if nargin >= 4 && ~isempty(n_segments_override)
        n_segments = n_segments_override;
    else
        n_segments = double(d.n_segments);
    end

    scale = (N / double(d.ref_N));    % keep the ring the same fraction of the canvas
    r_mid = double(d.r_mid) * scale;

    cx = (N+1)/2; cy = (N+1)/2;
    [Xg, Yg] = meshgrid(1:N, 1:N);
    dx = Xg - cx; dy = Yg - cy;
    theta = atan2(dy, dx);

    slot    = 2*pi / n_segments;
    seg_idx = round(theta / slot);
    theta_center = seg_idx * slot;
    tang_x = -sin(theta_center); tang_y = cos(theta_center);
    rad_x  =  cos(theta_center); rad_y  = sin(theta_center);

    seg_cx = cx + r_mid .* rad_x;
    seg_cy = cy + r_mid .* rad_y;
    ux = Xg - seg_cx; uy = Yg - seg_cy;
    u = (ux .* tang_x + uy .* tang_y) / scale;   % local coords, in template px
    v = (ux .* rad_x  + uy .* rad_y ) / scale;

    sampled = interp2(d.u_centers, d.v_centers, d.template, u, v, 'linear', d.meanL_ref);

    peak_dev  = max(d.template(:)) - d.meanL_ref;
    lum_scale = (contrast * meanL) / peak_dev;
    I_lum = meanL + (sampled - d.meanL_ref) * lum_scale;
end
