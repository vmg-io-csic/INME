function TCB = make_flash(N, n_frames, framerate, disc_mask, ...
                           onset_ms, dur_ms, meanL, amplitude, envelope_type)
% MAKE_FLASH  Generate a static luminance disc stimulus (flash onset/offset).
%
%   TCB = make_flash(N, n_frames, framerate, disc_mask,
%                    onset_ms, dur_ms, meanL, amplitude, envelope_type)
%
%   Creates a uniform-luminance disc that appears and disappears against
%   a uniform background.  Before onset and after offset, the entire
%   image is at meanL.  During the flash, pixels inside disc_mask are
%   set to meanL + amplitude.
%
% Inputs
%   N (int) image size [px] (square: N x N)
%   n_frames (int) total number of frames
%   framerate (double) frames per second [fps]
%   disc_mask (logical, N x N) binary mask from make_disc_mask
%   onset_ms (double) time when the disc appears [ms]
%   dur_ms (double) how long the disc stays visible [ms]
%   meanL (double) background luminance [a.u.]
%   amplitude (double) luminance increment inside the disc [a.u.]
%   envelope_type (char) temporal profile:
%                              'square'   — step onset/offset
%                              'gaussian' — smooth Gaussian envelope
%                                           (sigma = dur_ms / 4)
%
% Outputs
%   TCB (single, N x N x n_frames) stimulus tensor
%        Background pixels = meanL; disc pixels = meanL + amplitude * env(t).

    t_ms = (0:n_frames-1) / framerate * 1000;

    switch lower(envelope_type)
        case 'square'
            env = double(t_ms >= onset_ms & t_ms < (onset_ms + dur_ms));
        case 'gaussian'
            centre_ms = onset_ms + dur_ms/2;
            sigma_ms  = dur_ms / 4;
            env       = exp(-(t_ms - centre_ms).^2 / (2*sigma_ms^2));
        otherwise
            error('make_flash: unknown envelope_type ''%s''.', envelope_type);
    end

    base_frame = single(meanL * ones(N, N));
    TCB        = zeros(N, N, n_frames, 'single');
    for t = 1:n_frames
        frame = base_frame;
        frame(disc_mask) = meanL + amplitude * env(t);
        TCB(:,:,t) = frame;
    end
end
