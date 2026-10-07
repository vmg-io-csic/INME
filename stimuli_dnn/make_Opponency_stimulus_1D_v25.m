function mI = make_Opponency_stimulus_1D_v25(N, n_frames, framerate, dov, ...
            sigma_t_ms, speed_dps, bar_width, meanL, direction, inver)
% MAKE_OPPONENCY_STIMULUS_1D_V25  Drifting bar, framerate-correct.
%
%   Identical output to make_Opponency_stimulus_1D.m when framerate = 1000.
%   The original had
%
%       bar_center = mod(start_pos + speed_dps * t_vec(tt)/1000, dov);
%                                                          ^^^^  bug B1
%
%   with 1000 hardcoded instead of `framerate`, so the bar advanced
%   speed_dps/1000 deg per FRAME regardless of the sampling rate. It only
%   happened to be correct because Opponency was exported at 1000 fps; at
%   any other rate the bar's speed silently stopped being the nominal one.
%
% Inputs as make_Opponency_stimulus_1D.m
%
% Outputs
%   mI (double, n_frames x N)

    x_deg = (0:N-1) * dov / N;
    t_vec = (0:n_frames-1);

    t_center = n_frames / 2;
    sigma_fr = sigma_t_ms / 1000 * framerate;
    G_env = exp(-(t_vec - t_center).^2 / (2 * sigma_fr^2));
    mI = zeros(n_frames, N);

    start_pos = 0;
    for tt = 1:n_frames
        % Toroidal wrap so the bar stays on the canvas at any speed.
        bar_center = mod(start_pos + speed_dps * t_vec(tt)/framerate, dov);
        d = abs(x_deg - bar_center);
        d = min(d, dov - d);
        mI(tt, d <= (bar_width / 2)) = meanL * G_env(tt);
    end

    if direction == -1
        mI = fliplr(mI);
    end
    if inver == 1
        mI = meanL - mI;
    end
end
