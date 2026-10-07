%% run_TimeReversal_sweep_1D.m
%
% Sweeps Nm and Nw and reports the time-reversal symmetry-breaking metric
%  M = R[S] + R[Theta S],  which is zero for a
% detector that inverts its response exactly under time reversal.
%
% Both temporal poolings are reported from the same engine run: the flat
% flat mean over the presentation and the Tb-weighted mean
% used elsewhere in this work. The duration is set per stimulus to a whole
% number of temporal periods that is also comfortably longer than the
% filters, and the spatial phase is averaged.
%
% Inputs
%   Nm_list (1xM double) causal Tm durations to sweep [ms]
%   Nw_list (1xW double) Tw integration windows to sweep [ms]
%   dur_factor (double) duration as a multiple of Nm_max + Nw_max
%   n_phases (int) equispaced spatial phases to average
%   stim_subset (1xS int) 1 = LightEdges, 2 = UnevenSquareWave
%   sigma_c (double) width of the Gaussian readout filter c_s [px]
%
% Outputs (saved to results/timereversal_sweep_<stamp>.mat before printing):
%   M, M_sem, V_result, I_result, ratio_signed (MxWxSx2)
%   last index: 1 = flat mean, 2 = Tb weighted
%

clear; clc; close all;

cfg   = stINRF3_config('timereversal');
model = cfg.model;

% Full grid. The corner Nm = [300 500 750] by Nw = [500 750 1000] was
% probed separately for the uneven square wave, where M grows towards
% larger values of both parameters.
Nm_list = [30 50 75 100 125 200 300 500 750 1000];
Nw_list = [15 30 75 150 300 500];

dur_factor  = 2;
n_phases    = 15;   % 3 for the full grid, 15 or more for a focused one
stim_subset = [1 2];         % 1 = LightEdges, 2 = UnevenSquareWave
sigma_c     = 57.74;

N          = cfg.stim.N;
dov        = cfg.stim.dov;
framerate  = cfg.stim.framerate;
meanL      = cfg.stim.meanL;
contrast   = cfg.stim.contrast;

all_types = {'LightEdges', 'UnevenSquareWave'};
all_v     = [cfg.stim.v_LightEdges,      cfg.stim.v_UnevenSquareWave];
all_per   = [cfg.stim.period_LightEdges, cfg.stim.period_UnevenSquareWave];
all_fly   = [cfg.stim.ratio_fly_LightEdges, cfg.stim.ratio_fly_UnevenSquareWave];

stim_types = all_types(stim_subset);
v_list     = all_v(stim_subset);
per_list   = all_per(stim_subset);
fly_list   = all_fly(stim_subset);

out_dir = stinrf_results_dir(mfilename);
if isempty(out_dir), out_dir = pwd; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

n_m = numel(Nm_list);
n_w = numel(Nw_list);
n_s = numel(stim_types);

dur_min = dur_factor * (max(Nm_list) + max(Nw_list));
T_ms    = 1000 * per_list ./ v_list;
dur_ms  = ceil(dur_min ./ T_ms) .* T_ms;

j0        = round(N/2);
j_sensors = (j0-99):(j0+100);
phases    = 2*pi * (0:n_phases-1) / n_phases;

n_jobs = n_m * n_w * n_phases;
[im_g, iw_g, ip_g] = ndgrid(1:n_m, 1:n_w, 1:n_phases);
im_vec = im_g(:);  iw_vec = iw_g(:);  ip_vec = ip_g(:);

fprintf('\n==== TIME-REVERSAL SWEEP OVER (Nm, Nw) ====\n');
fprintf('%d x %d grid | %d equispaced phases | flat and Tb pooling | sigma_c = %.2f px\n', ...
        n_m, n_w, n_phases, sigma_c);
fprintf('filters need at least %g ms (Nm_max + Nw_max = %g ms)\n', ...
        dur_min, max(Nm_list) + max(Nw_list));
for is = 1:n_s
    fprintf('  %-18s period %.0f ms -> %.0f ms presentation = %.0f whole cycles\n', ...
            stim_types{is}, T_ms(is), dur_ms(is), dur_ms(is)/T_ms(is));
end

pool_names   = {'flat', 'Tb'};
M            = zeros(n_m, n_w, n_s, 2);
M_sem        = zeros(n_m, n_w, n_s, 2);
V_result     = zeros(n_m, n_w, n_s, 2);
I_result     = zeros(n_m, n_w, n_s, 2);
R_fwd_result = zeros(n_m, n_w, n_s, 2);
R_rev_result = zeros(n_m, n_w, n_s, 2);
ratio_signed = zeros(n_m, n_w, n_s, 2);

for is = 1:n_s
    stype = stim_types{is};
    v_k   = v_list(is);  per_k = per_list(is);
    n_frames = round(dur_ms(is) / 1000 * framerate);

    f_flat = zeros(1, n_jobs);   f_tbw = zeros(1, n_jobs);
    b_flat = zeros(1, n_jobs);   b_tbw = zeros(1, n_jobs);

    tn = (0:n_frames-1)' / n_frames;
    Tb = tn.^model.lpa .* (1 - tn).^model.lpb;
    Tb = Tb / (sum(Tb) + eps);

    parfor jj = 1:n_jobs
        m = model;
        m.Nm = Nm_list(im_vec(jj));
        m.Nw = Nw_list(iw_vec(jj));

        S_fwd = make_TimeReversal_stimulus_1D(stype, N, n_frames, framerate, ...
                    dov, phases(ip_vec(jj)), contrast, meanL, v_k, per_k);
        S_rev = flipud(S_fwd);

        [rf_fl, rf_tb] = both_poolings(S_fwd, dov, m, framerate, N, ...
                                       j_sensors, sigma_c, Tb);
        [rr_fl, rr_tb] = both_poolings(S_rev, dov, m, framerate, N, ...
                                       j_sensors, sigma_c, Tb);
        f_flat(jj) = rf_fl;   f_tbw(jj) = rf_tb;
        b_flat(jj) = rr_fl;   b_tbw(jj) = rr_tb;
    end

    for w = 1:2
        if w == 1
            Rf = reshape(f_flat, n_m, n_w, n_phases);
            Rr = reshape(b_flat, n_m, n_w, n_phases);
        else
            Rf = reshape(f_tbw,  n_m, n_w, n_phases);
            Rr = reshape(b_tbw,  n_m, n_w, n_phases);
        end
        Mp = Rf + Rr;
        M(:,:,is,w)     = mean(Mp, 3);
        M_sem(:,:,is,w) = std(Mp, 0, 3) / sqrt(n_phases);
        R_fwd_result(:,:,is,w) = mean(Rf, 3);
        R_rev_result(:,:,is,w) = mean(Rr, 3);
        V_result(:,:,is,w) = (R_fwd_result(:,:,is,w) - R_rev_result(:,:,is,w)) / 2;
        I_result(:,:,is,w) =  M(:,:,is,w) / 2;
        ratio_signed(:,:,is,w) = V_result(:,:,is,w) ./ sign_guard(I_result(:,:,is,w));
    end

    save(fullfile(out_dir, ['timereversal_sweep_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), ...
         'Nm_list','Nw_list','stim_types','dur_ms','dur_factor','n_phases','phases', ...
         'sigma_c','model','pool_names','M','M_sem','V_result','I_result', ...
         'R_fwd_result','R_rev_result','ratio_signed');

    fprintf('\n--- %s (v = %.2f deg/s, period = %.2f deg, %.0f ms, fly ratio %.2f) ---\n', ...
            stype, v_k, per_k, dur_ms(is), fly_list(is));

    for w = 1:2
        fprintf('\n  ======== pooling: %s ========\n', pool_names{w});
        print_grid('Symmetry-breaking metric  M = R[S] + R[Theta S]', ...
                   Nm_list, Nw_list, M(:,:,is,w), '%10.3f');
        print_grid('M in units of its SEM across phases', ...
                   Nm_list, Nw_list, M(:,:,is,w) ./ max(M_sem(:,:,is,w), eps), '%10.1f');
        print_grid('Signed ratio V / I', ...
                   Nm_list, Nw_list, ratio_signed(:,:,is,w), '%10.2f');

        sig = abs(M(:,:,is,w)) > 2 * M_sem(:,:,is,w);
        fprintf('  M differs from zero by more than 2 SEM in %d of %d combinations (%.0f%%)\n', ...
                sum(sig(:)), n_m*n_w, 100*sum(sig(:))/(n_m*n_w));
        pos = sig & (ratio_signed(:,:,is,w) > 0);
        fprintf('  of those, %d also have a correctly signed ratio\n', sum(pos(:)));
    end
end

fprintf('\nSaved to %s\n\n', out_dir);

function [r_flat, r_tb] = both_poolings(S, dov, m, framerate, N, j_sensors, sigma_c, Tb)
% BOTH_POOLINGS  Gaussian readout pooled flat and weighted by Tb.
    [~, F] = stINRF3_compute_1D(S, dov, m, framerate, N, j_sensors);
    n = size(F, 1);
    r = zeros(n, 1);
    for tt = 1:n
        r(tt) = stINRF3_readout_gaussian_1D(F(tt,:), j_sensors, sigma_c);
    end
    r_flat = mean(r);
    r_tb   = sum(Tb .* r);
end

function print_grid(header, Nm_list, Nw_list, G, fmt)
% PRINT_GRID  One metric as Nm rows by Nw columns.
    fprintf('\n  %s\n', header);
    fprintf('  %-10s', 'Nm \ Nw');
    for k = 1:numel(Nw_list), fprintf('%10s', sprintf('%d ms', Nw_list(k))); end
    fprintf('\n  %s\n', repmat('-', 1, 10 + 10*numel(Nw_list)));
    for a = 1:numel(Nm_list)
        fprintf('  %-10s', sprintf('%d ms', Nm_list(a)));
        fprintf(fmt, G(a,:));
        fprintf('\n');
    end
end

function s = sign_guard(x)
% SIGN_GUARD  Keeps the sign of x while avoiding division by zero.
    s = x;
    s(x == 0) = eps;
end
