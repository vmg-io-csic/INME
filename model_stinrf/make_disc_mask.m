function disc_mask = make_disc_mask(N, dov, radius_deg)
% MAKE_DISC_MASK  Generate a circular binary mask centred on the image.
%
% Inputs
%   N (int) image size [px] (square: N x N)
%   dov (double) field of view [degrees of visual angle]
%   radius_deg (double) disc radius [degrees of visual angle]
%
% Outputs
%   disc_mask (logical, N x N) true for pixels inside the disc

    dpi  = N / dov;
    r_px = round(radius_deg * dpi);
    [XX, YY] = meshgrid(1:N, 1:N);
    cx = (N+1)/2;  cy = (N+1)/2;
    disc_mask = ((XX-cx).^2 + (YY-cy).^2) <= r_px^2;
end
