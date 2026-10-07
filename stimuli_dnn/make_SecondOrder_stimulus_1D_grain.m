function mI = make_SecondOrder_stimulus_1D_grain(N, n_frames, framerate, dov, ...
        sigma_t_ms, phase_mov, phase_carrier, speed_dps, sfm_envelope_cpd, ...
        sfm_carrier_cpd, meanL, Type, Carrier_type, carrier_grain_px)
% MAKE_SECONDORDER_STIMULUS_1D_GRAIN  Second-order motion, coarse carrier.
%
%   Same as make_SecondOrder_stimulus_1D.m but the static binary carrier is
%   built in blocks of carrier_grain_px pixels.
%
%   WHY. The original carrier is 1 px binary noise, which RAFT, PWC-Net and
%   LiteFlowNet2 erase in their first x8 downsampling: the network never
%   sees the texture whose contrast the envelope modulates, so the stimulus
%   degenerates into a flat field before any motion computation happens.
%   A grain of 16 px (2x the x8 stride) survives it.
%
%   The carrier stays STATIC. There is no first-order motion energy in this
%   stimulus by construction, which is the point of the experiment: a
%   network that only tracks luminance features should report nothing.
%
% Inputs
%   carrier_grain_px (int) width of each carrier block [px]; 1 reproduces
%                          the original generator
%   (the rest as make_SecondOrder_stimulus_1D.m)
%
% Outputs
%   mI (double, n_frames x N)

    if nargin < 14 || isempty(carrier_grain_px), carrier_grain_px = 1; end

    x_deg = (0:N-1) * dov / N;
    t_vec = 0 : n_frames-1;

    if Carrier_type == 1
        n_blocks   = ceil(N / carrier_grain_px);
        block_vals = 2 * (rand(1, n_blocks) > 0.5) - 1;
        Carrier    = repelem(block_vals, carrier_grain_px);
        Carrier    = Carrier(1:N);
    else
        Carrier = sin(2*pi*sfm_carrier_cpd*x_deg + phase_carrier);
    end

    if Type == -1, direction = -1; else, direction = 1; end

    mI = zeros(n_frames, N);
    for tt = 1:n_frames
        Envelope = sin(2*pi*sfm_envelope_cpd * ...
                       (x_deg - direction*speed_dps*(t_vec(tt)/framerate)) ...
                       + phase_mov);
        mI(tt,:) = meanL * (1 + 0.5*(1 + Envelope) .* Carrier);
    end
end
