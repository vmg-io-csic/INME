%% run_Gratings.m
%
% Direction reported for a drifting grating, against the direction the
% stimulus specifies.
%
% One phase, not a random set. With a single phase there is no dispersion to
% report, so the figure carries no error bars and what it claims is that the
% model reports this direction, not that it reports it to a given precision.
% If a precision is wanted, sweep phases and report a circular SD.
%
% Inputs
%   grating_dirs (1xK double) directions to test [deg]
%   grating_sf (double) spatial frequency [cpd]
%   grating_speed_hz (double) temporal frequency [Hz]
%   half_region, readout_mode  readout
%
% Outputs:
%   results/Gratings/gratings_<stamp>.mat
%   export/Gratings/  via export_Gratings
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
grating_dirs     = 0:15:90;
grating_sf       = 0.25;
grating_speed_hz = 10;
phase            = 0;

%% READOUT
half_region  = 100;
readout_mode = 'gauss';

%% DERIVED
n_frames = round(duration_ms/1000*framerate);
nG       = numel(grating_dirs);
cfg      = stINRF3_config('gratings');  model = cfg.model;
f_cpx    = grating_sf / (N/dov);
sp       = grating_speed_hz / framerate;

fprintf('\n==== GRATINGS ====\n');
fprintf('%d directions | %.2f cpd | %g Hz | phase %g | readout %s\n', ...
        nG, grating_sf, grating_speed_hz, phase, readout_mode);

%% RUN
pred = zeros(nG, 1);
parfor k = 1:nG
    TCB = make_grating(N, n_frames, f_cpx, sp, grating_dirs(k), phase, ...
                       contrast, meanL);
    [O_H, O_V] = stINRF3_compute(TCB, dov, model, framerate, N);
    pred(k) = stinrf_readout_dir(O_H, O_V, readout_mode, half_region);
end

err = stINRF3_angular_error(pred, grating_dirs(:));
fprintf('\n  %-10s %8s %8s\n', 'Gratings', 'MAE', 'worst');
fprintf('  %-10s %8.2f %8.2f\n', '', mean(abs(err)), max(abs(err)));

%% SAVE
out_dir = stinrf_results_dir(mfilename);
save(fullfile(out_dir, ['gratings_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), ...
     'grating_dirs', 'grating_sf', 'grating_speed_hz', 'phase', 'pred', ...
     'err', 'N', 'dov', 'framerate', 'duration_ms', 'readout_mode', ...
     'half_region', 'model');
fprintf('Saved to %s\n', out_dir);

export_Gratings;
