%% run_Contrast_1D.m
%
% Runs the contrast saturation experiment
%
% Each condition is run once with the spatial phase fixed at 0, the same
% phase used for the optical-flow networks, so that both see identical
% stimuli. The model response is phase-invariant here (under 0.2% change
% across phases), so a single fixed phase loses nothing and is reproducible.
% The output is the mean model response per condition (signed, rightward
% positive).
%
% Inputs
%   N            (int)    number of spatial pixels
%   dov          (double) field of view [dva]
%   framerate    (int)    frame rate [fps]
%   duration_ms  (int)    total duration [ms]
%   speed_dps    (double) grating drift speed in degrees/second
%   meanL        (double) mean luminance [a.u.]
%   sfm_list     (1)      moving component SF per condition in cpd
%   sigma_t_list (double) Gaussian envelope width to sweep in ms
%   condition    (int)    condition index 1-5 (see CONDITION SELECTOR)
%   n_trials     (int)    repetitions per condition (1: phase is fixed)
%   rng_seed     (int)    random seed, 0 for unseeded
%   Contrast     (int)    Contrast of the grating
%

clear; close all;

cfg   = stINRF3_config('Cont');
model = cfg.model;

%% STIMULUS PARAMETERS
N            = 400;          % spatial samples [px]
dov          = 4.0;          % field of view [dva]
framerate    = 120;         % temporal sampling [fps] 
duration_ms  = 500;          % stimulus duration [ms]
speed_dps    = 8;           % grating drift speed [degrees/second]
meanL        = 100;          % mean luminance [a.u.]

% Spatial frequency [cpd]
sfm_list = [1]; 

Contrasts = logspace(-2,0,20);%%% 20 equidistant contrast in log scale

% Gaussian temporal envelope widths to sweep [ms]
sigma_t=100;


%% TRIAL / READOUT
n_trials = 1;               % one run, phase fixed at 0 (as for the DNNs)
rng_seed = 0;                % 0 = no seed (different each run)

%% READOUT
% Spatial readout filter c_s: central row of the 2D isotropic Gaussian,
% s_min = s_max.  sigma_c = L/sqrt(12) with L = 200 px the width of the
% sensor array it replaces.
sigma_c      = 57.74;        % [px]
make_videos  = false;        % true = also write the stimulus videos

%% DERIVED
n_frames   = round(duration_ms / 1000 * framerate);
n_contrast    = numel(Contrasts);
j0         = round(N/2);
j_sensors  = (j0-99):(j0+100);  % 200 sensors centred on the stimulus

if rng_seed > 0, rng(rng_seed); end

% Pre-calculate physical units for axes labeling
x_deg = (0:N-1) * dov / N - dov/2;
t_sec = (0:n_frames-1) / framerate;

fprintf('\n==== RESULTS CONTRAST SATURATION (1D) ====\n');
fprintf('Readout: 200 sensors\n');

%% RUN
fprintf('  %-12s %-16s %-12s\n', 'Contrast', 'Response', 'Trials');
fprintf('  %s\n', repmat('-', 1, 44));

for si = 1:n_contrast
    n_correct = 0;
    for tr = 1:n_trials
        phase_mov  = 0;              % fixed phase, as for the DNNs

        mI = make_Contrast_stimulus_1D(N, n_frames, framerate, dov, ...
            sigma_t, phase_mov, ...
            speed_dps, sfm_list, meanL);


 mI=meanL*(1+Contrasts(si)*2*(mI./abs(max(max(mI)))-0.5));

 figure(1)       
 imagesc(x_deg, t_sec, mI, [0 200])
 colormap(gray)
 axis square
 xlabel('Space (degrees)')
 ylabel('Time (seconds)')
 pause(0.1)


        O_sensors = stINRF3_compute_1D(mI, dov, model, framerate, N, j_sensors);

       %figure(2);hist(O_sensors);pause(0.01)
       Response(tr) = stINRF3_readout_gaussian_1D(O_sensors, j_sensors, sigma_c);

    end

    Mean_response(si) = mean(Response);

    fprintf('  %-12s %-16s\n', ...
        sprintf('%.3f', Contrasts(si)), ...
        sprintf('%d    %d', mean(Response), n_trials));
end

fprintf('\n');

figure(2)
loglog(100*Contrasts,Mean_response,'bo-')
ax = gca; 
ax.XTick = [1, 3, 10, 30, 100]; 
ax.XTickLabel = {'1', '3', '10', '30', '100'}; 
ax.YTick = [0.001, 0.01, 0.1, 1, 2, 5]; 
ax.YTickLabel = {'10^{-3}', '0.01', '0.1', '1', '2', '5'};
ylim([0.001 10]); 
ylabel('MS-INRF response (a.u.)')
xlabel('% Contrast')
axis square
grid on

%% SAVE  (small variables only, so the plotting script can reuse them)
res_dir = stinrf_results_dir(mfilename);
if isempty(res_dir), res_dir = pwd; end
if ~exist(res_dir,'dir'), mkdir(res_dir); end
wv = whos; small = {wv([wv.bytes] < 2e5).name};
save(fullfile(res_dir, ['Contrast_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), small{:});
fprintf('Saved to %s\n', res_dir);

export_Contrast;
