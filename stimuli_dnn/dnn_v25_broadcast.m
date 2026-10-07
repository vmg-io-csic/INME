function tensor = dnn_v25_broadcast(mI, N_rows)
% DNN_V25_BROADCAST  1D stimulus (frames x cols) -> (rows x cols x frames).
    n_frames = size(mI, 1);
    N_cols   = size(mI, 2);
    tensor = zeros(N_rows, N_cols, n_frames, 'single');
    for t = 1:n_frames
        tensor(:, :, t) = single(repmat(mI(t, :), [N_rows, 1]));
    end
end
