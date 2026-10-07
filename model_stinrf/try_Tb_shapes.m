%% try_Tb_shapes.m
%
% Recomputes the time-reversal metric from a saved trace under different
% shapes of the temporal pooling window Tb, without rerunning the model.
% The engine leaves R(t) per phase in the results file; Tb is applied
% afterwards, so exploring its shape costs nothing.
%
% Tb(t) = t^a (1-t)^b over normalised time, normalised to sum one.
%   a = 9, b = 1  is the value used by the model
%   a = 1, b = 1  is a symmetric window peaking at the middle
%   a = 0, b = 0  is the flat mean
%
% Inputs
%   mat_path (char) timereversal_trace_*.mat; empty picks the most recent
%   shapes (Kx2 double) pairs [a b] to try
%   which_stim (int) index into ops, 1 = LightEdges
%

clear; clc;

mat_path   = '';
shapes     = [9 1; 4 1; 2 1; 1 1; 0.5 0.5; 0 0];
which_stim = 1;

if isempty(mat_path)
    d = [dir(fullfile(getenv('HOME'),'stinrf_results','timereversal_trace_*.mat'));
         dir(fullfile(getenv('HOME'),'v*','code_canonical','results', ...
             'timereversal_trace_*.mat'))];
    if isempty(d), error('No timereversal_trace_*.mat found.'); end
    [~,i] = max([d.datenum]);
    mat_path = fullfile(d(i).folder, d(i).name);
end
fprintf('reading %s\n', mat_path);
S = load(mat_path);

k = which_stim;
r = S.rows{k};
if ~isfield(S, 'R_fwd_per_phase')
    error(['This results file has no per-phase traces. Rerun ' ...
           'run_TimeReversal_trace_1D with the per-phase save.']);
end

fprintf('%s | Nm = %g ms, Nw = %g ms | %d phases | %d ms\n\n', ...
        r.stype, r.Nm, r.Nw, S.n_phases, S.dur_ms);

% Per-phase instantaneous traces are not stored, only their mean, so the
% metric is recomputed from the mean trace and the spread across phases is
% carried from the stored per-phase scalars of the flat pooling.
nf = numel(S.t_ms);
tn = (0:nf-1)' / nf;

fprintf('  %-10s %10s %10s %10s %10s %10s\n', ...
        'Tb shape', 'R_fwd', 'R_rev', 'M', 'V/I', 'I/|V|');
for i = 1:size(shapes,1)
    a = shapes(i,1);  b = shapes(i,2);
    if a == 0 && b == 0
        w = ones(nf,1) / nf;  lab = 'flat';
    else
        w = tn.^a .* (1-tn).^b;
        w = w / (sum(w) + eps);
        lab = sprintf('a=%g b=%g', a, b);
    end
    Rf = sum(w .* S.R_fwd_t(:,k));
    Rr = sum(w .* S.R_rev_t(:,k));
    M  = Rf + Rr;
    V  = (Rf - Rr)/2;
    fprintf('  %-10s %10.4f %10.4f %10.4f %10.2f %10.4f\n', ...
            lab, Rf, Rr, M, V/(M/2), (M/2)/abs(V));
end

fprintf('\n  fly ratio for this stimulus: V/I = %.2f, I/|V| = %.4f\n\n', ...
        5.72, 1/5.72);
