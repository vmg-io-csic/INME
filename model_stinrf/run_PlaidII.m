%% run_PlaidII.m
%
% Direction reported for a Type-II plaid, against the intersection-of-
% constraints direction. Each configuration gives its two components a
% different speed and orientation, so the IOC direction falls outside the
% angle between them, which is what makes these the discriminating case.
%
% One phase, not a random set; see the note in run_Gratings.m.
%
% Inputs
%   plaidII_sf (double) component spatial frequency [cpd]
%   plaidII_configs (Kx4 double) speed1, speed2, theta1, theta2 per
%       configuration. The set below is symmetric in even 10 deg steps from
%       -30 to +30. Two earlier entries at +/-25.95 were dropped: they sat
%       four degrees from the ends, which forced their stimulus icons to half
%       the size of the ones in the other two figures, and they sampled
%       nothing the 10 deg grid does not already cover
%
% Outputs:
%   results/PlaidII/plaidII_<stamp>.mat
%   export/PlaidII/  via export_PlaidII
%

clear; clc; close all;

%% STIMULUS GEOMETRY
N           = 300;
dov         = 300*(1.5/128);
framerate   = 200;
duration_ms = 2000;
contrast    = 1;
meanL       = 100;

%% EXPERIMENT
plaidII_sf = 2.5;
plaidII_configs = [ ...
    0.6046, 1.2045,   42,   22;    % IOC -30.00
    0.6046, 1.2045,   52,   32;    % IOC -20.00
    0.6661, 1.0931,   54,   34;    % IOC -10.00
    0.7185, 1.0134,   54,   34;    % IOC   0.00
    0.6661, 1.0931,  -54,  -34;    % IOC  10.00
    0.6046, 1.2045,  -52,  -32;    % IOC  20.00
    0.6046, 1.2045,  -42,  -22];   % IOC  30.00
phase = 0;

%% READOUT
half_region  = 100;
readout_mode = 'gauss';

%% DERIVED
n_frames = round(duration_ms/1000*framerate);
nQ       = size(plaidII_configs, 1);
cfg      = stINRF3_config('plaidII');  model = cfg.model;
f_cpx    = plaidII_sf / (N/dov);
tf_scale = plaidII_sf / 0.21333 / framerate;

ioc = zeros(nQ, 1);
for c = 1:nQ
    ioc(c) = plaidII_ioc(plaidII_configs(c,1), plaidII_configs(c,2), ...
                         plaidII_configs(c,3), plaidII_configs(c,4));
end

fprintf('\n==== TYPE-II PLAIDS ====\n');
fprintf('%d configurations | %.2f cpd | phase %g | readout %s\n', ...
        nQ, plaidII_sf, phase, readout_mode);

%% RUN
sp1 = plaidII_configs(:,1) * tf_scale;
sp2 = plaidII_configs(:,2) * tf_scale;
th1 = plaidII_configs(:,3);
th2 = plaidII_configs(:,4);

pred = zeros(nQ, 1);
parfor k = 1:nQ
    TCB = make_plaidII(N, n_frames, f_cpx, sp1(k), sp2(k), th1(k), th2(k), ...
                       phase, phase, contrast, meanL);
    [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
    pred(k) = stinrf_readout_dir(O_H, O_V, readout_mode, half_region);
end

err = stINRF3_angular_error(pred, ioc);
fprintf('\n  %-10s %8s %8s\n', 'Plaid II', 'MAE', 'worst');
fprintf('  %-10s %8.2f %8.2f\n', '', mean(abs(err)), max(abs(err)));

%% SAVE
out_dir = stinrf_results_dir(mfilename);
save(fullfile(out_dir, ['plaidII_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), ...
     'plaidII_sf', 'plaidII_configs', 'ioc', 'phase', 'pred', 'err', ...
     'N', 'dov', 'framerate', 'duration_ms', 'readout_mode', ...
     'half_region', 'model');
fprintf('Saved to %s\n', out_dir);

export_PlaidII;
