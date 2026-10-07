function TCB = make_grating_ref(N, n_frames, f_cpx, speed, meanL)
% MAKE_GRATING_REF  Generate a rightward drifting grating for magnitude reference.
%
%   TCB = make_grating_ref(N, n_frames, f_cpx, speed, meanL)
%
%   Creates a horizontal (0°) drifting grating with full contrast and
%   zero initial phase.  Used to compare the model's response magnitude
%   on the flash disc ring against its response to real drifting motion
%   with the same spatial/temporal parameters.
%
% Inputs
%   N (int) image size [px] (square: N x N)
%   n_frames (int) number of temporal frames
%   f_cpx (double) spatial frequency [cycles/pixel]
%   speed (double) temporal frequency [cycles/frame]
%   meanL (double) mean luminance [a.u.]
%
% Outputs
%   TCB (single, N x N x n_frames) stimulus tensor (horizontal grating)

    coords = linspace(-N/2, N/2, N);
    [X, ~] = meshgrid(coords, coords);

    TCB = zeros(N, N, n_frames, 'single');
    for t = 1:n_frames
        phase_t    = 2*pi * speed * (t-1);
        TCB(:,:,t) = single(meanL * (1 + sin(2*pi * f_cpx * X + phase_t)));
    end
end
