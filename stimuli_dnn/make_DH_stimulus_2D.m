function TCB = make_DH_stimulus_2D(N, N_rows, n_frames, framerate, dov, ...
                                    cond_idx, sigma_t_ms, phase_mov, phase_stat, ...
                                    speed_dps, sfm_list, sfs_list, meanL)
% MAKE_DH_STIMULUS_2D  Generate a D&H stimulus as a 2D image (replicated rows).
%
%   TCB = make_DH_stimulus_2D(N, N_rows, n_frames, framerate, dov,
%                              cond_idx, sigma_t_ms, phase_mov, phase_stat,
%                              speed_dps, sfm_list, sfs_list, meanL)
%
%   Wraps make_DH_stimulus_1D: generates the 1D row stimulus and broadcasts
%   it to all N_rows rows, producing a stimulus that is uniform in the
%   vertical direction.  This allows using the full 2D model engine
%   (stINRF3_compute) on D&H stimuli.
%
% Inputs
%   N, N_rows, n_frames, framerate, dov, cond_idx, sigma_t_ms,
%   phase_mov, phase_stat, speed_dps, sfm_list, sfs_list, meanL
%       — Same as make_DH_stimulus_1D, plus:
%   N_rows (int) number of image rows (height). Typically 50
%
% Outputs
%   TCB (single, N_rows x N x n_frames) stimulus tensor

    mI  = make_DH_stimulus_1D(N, n_frames, framerate, dov, cond_idx, ...
                               sigma_t_ms, phase_mov, phase_stat, ...
                               speed_dps, sfm_list, sfs_list, meanL);

    TCB = zeros(N_rows, N, n_frames, 'single');
    for tt = 1:n_frames
        TCB(:,:,tt) = repmat(mI(tt,:), [N_rows, 1]);
    end
end
