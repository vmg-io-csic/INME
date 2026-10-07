function stim_1D = make_TimeReversal_stimulus_1D(stim_type, N, n_frames, ...
                                                  framerate, dov, phase0_rad, ...
                                                  contrast, meanL, ...
                                                  v_dps, period_deg)
% MAKE_TIMEREVERSAL_STIMULUS_1D  Time-reversal symmetry-breaking stimulus.
%
%   stim_1D = make_TimeReversal_stimulus_1D(stim_type, N, n_frames, ...
%                 framerate, dov, phase0_rad, contrast, meanL, ...
%                 v_dps, period_deg)
%
%   Two periodic 1D profiles, each with a different degree of temporal
%   asymmetry, drifting at constant velocity v_dps:
%     'LightEdges'       — single sharp transition per period:
%                           S(x,t) = +1 if (x mod period) <= (v*t mod period),
%                           -1 otherwise.
%     'UnevenSquareWave' — uneven square wave, relative bar widths
%                           5:5:25:25 in units of period/60.
%
%   Reversing the frame order of the returned stimulus (time reversal, NOT
%   flipping the sign of v_dps) probes whether a motion detector is
%   perfectly antisymmetric under time reversal. A biological or model
%   detector that is NOT perfectly antisymmetric produces a nonzero
%   symmetric component when forward and reversed responses are combined
%   — see run_TimeReversal_1D.m for the V/I/ratio bookkeeping.
%
%   Defaults if v_dps / period_deg are omitted (or empty), matching the
%   paper's own stimuli:
%     LightEdges:       v_dps = 9.38,  period_deg = 1.17
%     UnevenSquareWave: v_dps = 1.17,  period_deg = 2.34
%
% Inputs
%   stim_type (char) 'LightEdges' or 'UnevenSquareWave'
%   N (int) number of spatial samples (pixels along x)
%   n_frames (int) number of temporal frames
%   framerate (double) frames per second [fps]
%   dov (double) field of view [degrees]
%   phase0_rad (double) initial spatial phase [radians]
%   contrast (double) Michelson-like contrast multiplier on the ±1 profile
%   meanL (double) mean luminance [a.u.]
%   v_dps (double, optional) drift speed [degrees/second]
%   period_deg (double, optional) spatial period [degrees]
%
% Outputs
%   stim_1D (double, n_frames x N) 1D stimulus row repeated across time

    if nargin < 9 || isempty(v_dps)
        switch stim_type
            case 'LightEdges',       v_dps = 9.38;
            case 'UnevenSquareWave', v_dps = 1.17;
            otherwise, error('make_TimeReversal_stimulus_1D: unknown stim_type ''%s''.', stim_type);
        end
    end
    if nargin < 10 || isempty(period_deg)
        switch stim_type
            case 'LightEdges',       period_deg = 1.17;
            case 'UnevenSquareWave', period_deg = 2.34;
        end
    end

    dt    = 1 / framerate;
    x_deg = linspace(-dov/2, dov/2, N);
    t_s   = (0:n_frames-1) * dt;

    c = zeros(n_frames, N);

    switch stim_type

        case 'LightEdges'
            % Single sharp transition per period, advancing with v_dps.
            x_shift = phase0_rad * period_deg / (2*pi);
            for tt = 1:n_frames
                pos_in_period = mod(x_deg + x_shift, period_deg) / period_deg;
                edge_position = mod(v_dps * t_s(tt), period_deg) / period_deg;
                row = -ones(1, N);
                row(pos_in_period <= edge_position) = 1;
                c(tt,:) = row;
            end

        case 'UnevenSquareWave'
            % Relative bar widths 5:5:25:25 (boundaries at 5/60, 10/60,
            % 35/60 of period_deg).
            for tt = 1:n_frames
                xs = mod((x_deg - v_dps*t_s(tt) + ...
                         phase0_rad*period_deg/(2*pi)) / period_deg, 1);
                row = -ones(1, N);
                row(xs >= 5/60  & xs < 10/60) = 1;
                row(xs >= 35/60 & xs < 1)     = 1;
                c(tt,:) = row;
            end

        otherwise
            error('make_TimeReversal_stimulus_1D: unknown stim_type ''%s'' (use ''LightEdges'' or ''UnevenSquareWave'').', stim_type);
    end

    stim_1D = meanL * (1 + contrast * c);
end
