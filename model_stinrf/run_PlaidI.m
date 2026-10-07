%% run_PlaidI.m
%
% Direction reported for a Type-I plaid, against the direction the pattern
% moves in. The two components sit at plus and minus plaidI_offset around it.
%
% One phase, not a random set; see the note in run_Gratings.m.
%
% Inputs
%   plaidI_dirs (1xK double) pattern directions to test [deg]
%   plaidI_sf (double) component spatial frequency [cpd]
%   plaidI_speed_hz (double) component temporal frequency [Hz]
%   plaidI_offset (double) half angle between the components [deg]
%
% Outputs:
%   results/PlaidI/plaidI_<stamp>.mat
%   export/PlaidI/  via export_PlaidI
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
plaidI_dirs     = 0:15:90;
plaidI_sf       = 1.0;
plaidI_speed_hz = 8;
plaidI_offset   = 45;
phase           = 0;

%% READOUT
half_region  = 100;
readout_mode = 'gauss';

%% DERIVED
n_frames = round(duration_ms/1000*framerate);
nP       = numel(plaidI_dirs);
cfg      = stINRF3_config('plaidI');  model = cfg.model;
f_cpx    = plaidI_sf / (N/dov);
sp       = plaidI_speed_hz / framerate;

fprintf('\n==== TYPE-I PLAIDS ====\n');
fprintf('%d directions | %.2f cpd | %g Hz | components at +/-%g deg | phase %g\n', ...
        nP, plaidI_sf, plaidI_speed_hz, plaidI_offset, phase);

%% RUN
pred = zeros(nP, 1);
parfor k = 1:nP
    TCB = make_plaidI(N, n_frames, f_cpx, sp, ...
                      plaidI_dirs(k) + plaidI_offset, ...
                      plaidI_dirs(k) - plaidI_offset, ...
                      phase, phase, contrast, meanL);
    [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
    pred(k) = stinrf_readout_dir(O_H, O_V, readout_mode, half_region);
end

err = stINRF3_angular_error(pred, plaidI_dirs(:));
fprintf('\n  %-10s %8s %8s\n', 'Plaid I', 'MAE', 'worst');
fprintf('  %-10s %8.2f %8.2f\n', '', mean(abs(err)), max(abs(err)));

%% SAVE
out_dir = stinrf_results_dir(mfilename);
save(fullfile(out_dir, ['plaidI_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), ...
     'plaidI_dirs', 'plaidI_sf', 'plaidI_speed_hz', 'plaidI_offset', ...
     'phase', 'pred', 'err', 'N', 'dov', 'framerate', 'duration_ms', ...
     'readout_mode', 'half_region', 'model');
fprintf('Saved to %s\n', out_dir);

export_PlaidI;
