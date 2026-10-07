%% plot_lightedges_trace.m
%
% Response of the model to the light-edge stimulus played forward and in
% reverse, built from the trace
% results file. Nothing is recomputed.
%
% Top row: the stimulus in space and time, framed in the colour of the
% direction it moves.
% Middle: the instantaneous readout averaged over spatial phases, with a
% shaded band for the standard deviation across phases.
% Bottom: the mean over the presentation, one dot per phase, with the
% standard deviation.
%
% Both temporal poolings are shown side by side: the flat mean over the
% presentation, and the Tb-weighted mean, which is the
% response the engine itself returns.
%
% Inputs
%   mat_path (char) timereversal_trace_*.mat to read; empty picks the most
%       recent one found under the home directory
%   which_stim (int) index into ops, 1 = LightEdges
%   space_win (double) half-width of the space axis in the panels [deg]
%   time_win (double) span of the space-time panels [ms]
%
% Outputs:
%   Fig_LightEdgesTrace.png at 600 dpi and .pdf, next to the results
%

clear; clc; close all;

mat_path   = '';
which_stim = 1;
space_win  = 3.6;
time_win   = 250;

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
out_dir = fileparts(mat_path);

k = which_stim;
r = S.rows{k};
fprintf('%s | Nm = %g ms, Nw = %g ms | %d %s phases | %d ms\n', ...
        r.stype, r.Nm, r.Nw, S.n_phases, S.phase_mode, S.dur_ms);
fprintf('\n  %-6s %12s %12s %12s %10s\n', 'pool','R_fwd','R_rev','M','A = I/|V|');
for w = 1:2
    if w == 1
        Rf = r.R_fwd_flat;  Rr = r.R_rev_flat;  nm = 'flat';
    else
        Rf = r.R_fwd_tb;    Rr = r.R_rev_tb;    nm = 'Tb';
    end
    M = Rf + Rr;  V = (Rf - Rr)/2;
    fprintf('  %-6s %12.4f %12.4f %12.4f %10.4f\n', nm, Rf, Rr, M, (M/2)/abs(V));
end

c = stinrf_fig_style();
cfg = stINRF3_config('timereversal');
N   = cfg.stim.N;  dov = cfg.stim.dov;  fr = cfg.stim.framerate;
nf  = min(size(S.stim_fwd{k},1), round(time_win/1000*fr));
xd  = (0:N-1)*dov/N;
keep = xd <= space_win;
tp  = (0:nf-1)/fr*1000;

fig = figure('Color','w','Units','inches','Position',[1 1 10.0 6.4], ...
             'PaperPositionMode','auto');

% ---- stimulus panels
for q = 1:2
    ax = axes('Units','normalized','Position',[0.06+(q-1)*0.175 0.60 0.145 0.33]);
    if q == 1
        img = S.stim_fwd{k};  col = c.right;  ttl = 'Time forward';
    else
        img = S.stim_rev{k};  col = c.left;   ttl = 'Time reversed';
    end
    imagesc(xd(keep), tp, img(1:nf, keep));
    colormap(gray); set(gca,'YDir','reverse');
    xlabel('Space (deg)','FontName','Arial','FontSize',9);
    if q == 1, ylabel('Time (ms)','FontName','Arial','FontSize',9); end
    stinrf_fig_style(ax);
    xl = xlim(ax);  yl = ylim(ax);
    rectangle(ax,'Position',[xl(1) yl(1) diff(xl) diff(yl)], ...
              'EdgeColor',col,'LineWidth',2.4);
    title(ax, ttl, 'FontWeight','normal','FontName','Arial','FontSize',9,'Color',col);
end

% ---- instantaneous readout, mean over phases with an SD band
t_ms = S.t_ms;
for w = 1:2
    ax = axes('Units','normalized','Position',[0.46+(w-1)*0.28 0.60 0.22 0.33]);
    hold(ax,'on');
    if w == 1
        yf = S.R_fwd_t(:,k);      yr = S.R_rev_t(:,k);
        sf = S.R_fwd_sem(:,k)*sqrt(S.n_phases);
        sr = S.R_rev_sem(:,k)*sqrt(S.n_phases);
        ttl = 'Instantaneous response';
    else
        yf = S.R_fwd_tb_t(:,k);   yr = S.R_rev_tb_t(:,k);
        sf = S.Tb .* S.R_fwd_sem(:,k)*sqrt(S.n_phases);
        sr = S.Tb .* S.R_rev_sem(:,k)*sqrt(S.n_phases);
        ttl = 'Weighted by T_b';
    end
    fill(ax, [t_ms; flipud(t_ms)], [yf+sf; flipud(yf-sf)], c.right, ...
         'FaceAlpha', 0.22, 'EdgeColor','none');
    fill(ax, [t_ms; flipud(t_ms)], [yr+sr; flipud(yr-sr)], c.left, ...
         'FaceAlpha', 0.22, 'EdgeColor','none');
    plot(ax, t_ms, yf, '-', 'Color', c.right, 'LineWidth', 1.2);
    plot(ax, t_ms, yr, '-', 'Color', c.left,  'LineWidth', 1.2);
    yline(ax, 0, 'k-', 'LineWidth', 0.8);
    xlabel('Time (ms)','FontName','Arial','FontSize',9);
    ylabel('R(t) (a.u.)','FontName','Arial','FontSize',9);
    stinrf_fig_style(ax, ttl);
    if w == 1
        yl = ylim(ax);
        text(ax, t_ms(end), yl(2), 'Forward ', 'Color', c.right, ...
             'FontName','Arial','FontSize',9,'FontWeight','bold', ...
             'HorizontalAlignment','right','VerticalAlignment','top');
        text(ax, t_ms(end), yl(1), 'Reversed ', 'Color', c.left, ...
             'FontName','Arial','FontSize',9,'FontWeight','bold', ...
             'HorizontalAlignment','right','VerticalAlignment','bottom');
    end
end

% ---- mean over the presentation, one dot per phase
for w = 1:2
    ax = axes('Units','normalized','Position',[0.46+(w-1)*0.28 0.11 0.22 0.33]);
    hold(ax,'on');
    if w == 1
        Rf = S.R_fwd_per_phase(:,k);  Rr = S.R_rev_per_phase(:,k);  ttl = 'Flat mean';
    else
        Rf = S.R_fwd_tb_per_phase(:,k);  Rr = S.R_rev_tb_per_phase(:,k);
        ttl = 'T_b weighted mean';
    end
    plot_runs_bar(ax, 1, Rf, c.right);
    plot_runs_bar(ax, 2, Rr, c.left);
    yline(ax, 0, 'k-', 'LineWidth', 0.8);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Forward','Reversed'});
    xlim(ax, [0.4 2.6]);
    ylabel('Mean response (a.u.)','FontName','Arial','FontSize',9);
    stinrf_fig_style(ax, ttl);
end

exportgraphics(fig, fullfile(out_dir,'Fig_LightEdgesTrace.png'),'Resolution',600);
exportgraphics(fig, fullfile(out_dir,'Fig_LightEdgesTrace.pdf'),'ContentType','vector');
fprintf('\nSaved Fig_LightEdgesTrace in %s\n\n', out_dir);
