function out = dnn_v25_crop(TCB, N_rows)
% DNN_V25_CROP  Central row crop of a square tensor to N_rows.
%
%   The 2D generators (make_grating, make_plaidI, make_plaidII) build a
%   square N x N canvas. The networks' canvases are rectangular, so the
%   stimulus is generated at N = N_cols and the rows are cropped centrally.
%   Because the generators lay x and y on the same pixel grid, cropping
%   rows changes the field of view in y without touching px_per_deg, the
%   spatial frequency or the speed.
    N = size(TCB, 1);
    if N_rows >= N, out = TCB; return; end
    r0 = floor((N - N_rows)/2) + 1;
    out = TCB(r0:r0+N_rows-1, :, :);
end
