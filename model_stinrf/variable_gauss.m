function out = variable_gauss(img, sigma_map, n_levels)
% VARIABLE_GAUSS  Spatially-variable Gaussian blur, discretised into
%   n_levels sigma bins for speed (one imgaussfilt call per bin instead
%   of one Gaussian per pixel).
%
%   out = variable_gauss(img, sigma_map, n_levels)
%
% Inputs
%   img (double, H x W) field to blur
%   sigma_map (double, H x W) per-pixel Gaussian sigma [px]
%   n_levels (int, optional) number of discrete sigma bins (default 12)
%
% Outputs
%   out (double, H x W) blurred field

    if nargin < 3 || isempty(n_levels), n_levels = 12; end

    s_min = min(sigma_map(:));
    s_max = max(sigma_map(:));

    if (s_max - s_min) < 1e-6
        out = imgaussfilt(img, s_min);
        return;
    end

    levels = linspace(s_min, s_max, n_levels);
    out = zeros(size(img));

    % Nearest discrete level per pixel
    [~, idx_map] = min(abs(reshape(sigma_map, [], 1) - levels), [], 2);
    idx_map = reshape(idx_map, size(sigma_map));

    for L = 1:n_levels
        m = (idx_map == L);
        if ~any(m(:)), continue; end
        blurred = imgaussfilt(img, max(levels(L), 1e-3));
        out(m) = blurred(m);
    end
end
