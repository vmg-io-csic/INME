%% plot_timereversal_summary.m
%
% Summary figure for the time-reversal experiments, in the layout of the
% response-trace figure: the two stimuli in space and time on the left, and
% on the right how the symmetry-breaking metric and the ratio depend on how
% the response is pooled over time and on how long the stimulus lasts.
%
% The point of the figure: M = R[S] + R[Theta S] is a small difference
% between two large numbers, so it is sensitive to the pooling and to the
% duration. Reading a single number without stating both is meaningless.
%
% Inputs
%   mat_path (char) timereversal_duration_*.mat to read; empty picks the
%       most recent one found under the home directory
%   space_win (double) half-width of the stimulus space axis [deg]
%   time_win (double) span of the space-time panels [ms]
%
% Outputs:
%   Fig_TimeReversalSummary.png at 600 dpi and .pdf, next to the results
%

clear; clc; close all;

mat_path  = '';
space_win = 3.6;
time_win  = 250;

if isempty(mat_path)
    d = [dir(fullfile(getenv('HOME'),'stinrf_results','timereversal_duration_*.mat'));
         dir(fullfile(getenv('HOME'),'v*','code_canonical','results', ...
             'timereversal_duration_*.mat'))];
    if isempty(d), error('No timereversal_duration_*.mat found.'); end
    [~,i] = max([d.datenum]);
    mat_path = fullfile(d(i).folder, d(i).name);
end
fprintf('reading %s\n', mat_path);
S = load(mat_path);
out_dir = fileparts(mat_path);

c = stinrf_fig_style();
cfg = stINRF3_config('timereversal');
N   = cfg.stim.N;  dov = cfg.stim.dov;  fr = cfg.stim.framerate;
v_of   = struct('LightEdges', cfg.stim.v_LightEdges, ...
                'UnevenSquareWave', cfg.stim.v_UnevenSquareWave);
per_of = struct('LightEdges', cfg.stim.period_LightEdges, ...
                'UnevenSquareWave', cfg.stim.period_UnevenSquareWave);

n_ops = numel(S.ops);
nf    = round(time_win/1000*fr);
x     = (0:N-1)*dov/N;
keep  = x <= space_win;
t     = (0:nf-1)/fr*1000;

fprintf('\n  %-18s %-6s %8s %10s %10s %8s\n', ...
        'stimulus','pool','dur','M','SEM(M)','V/I');
for k = 1:n_ops
    for w = 1:2
        for id = 1:numel(S.dur_list)
            fprintf('  %-18s %-6s %6d ms %10.4f %10.4f %8.2f\n', ...
                    S.ops(k).stype, S.wnames{w}, S.dur_list(id), ...
                    S.M(k,id,w), S.M_sem(k,id,w), S.ratio(k,id,w));
        end
    end
end

fig = figure('Color','w','Units','inches','Position',[1 1 10.5 5.6], ...
             'PaperPositionMode','auto');

for k = 1:n_ops
    stype = S.ops(k).stype;
    yb = 0.60 - (k-1)*0.47;

    Sf = make_TimeReversal_stimulus_1D(stype, N, nf, fr, dov, 0, ...
             cfg.stim.contrast, cfg.stim.meanL, v_of.(stype), per_of.(stype));
    Sr = flipud(Sf);

    for q = 1:2
        ax = axes('Units','normalized','Position',[0.045+(q-1)*0.135 yb 0.11 0.31]);
        if q == 1, img = Sf; col = c.right; ttl = 'Time forward';
        else,      img = Sr; col = c.left;  ttl = 'Time reversed'; end
        imagesc(x(keep), t, img(:,keep));
        colormap(gray); set(gca,'YDir','reverse');
        xlabel('Space (deg)','FontName','Arial','FontSize',9);
        if q == 1, ylabel('Time (ms)','FontName','Arial','FontSize',9); end
        stinrf_fig_style(ax);
        xl = xlim(ax);  yl = ylim(ax);
        rectangle(ax,'Position',[xl(1) yl(1) diff(xl) diff(yl)], ...
                  'EdgeColor',col,'LineWidth',2.4);
        title(ax, ttl, 'FontWeight','normal','FontName','Arial', ...
              'FontSize',9,'Color',col);
    end

    ax = axes('Units','normalized','Position',[0.375 yb 0.24 0.31]); hold on;
    errorbar(S.dur_list, S.M(k,:,1), S.M_sem(k,:,1), '-o', 'Color', c.right, ...
             'LineWidth',1.6,'MarkerFaceColor','w','MarkerSize',6,'CapSize',5);
    errorbar(S.dur_list, S.M(k,:,2), S.M_sem(k,:,2), '-s', 'Color', c.ref, ...
             'LineWidth',1.6,'MarkerFaceColor','w','MarkerSize',6,'CapSize',5);
    yline(0,'k--','LineWidth',1);
    set(gca,'XTick',S.dur_list);
    xlabel('Duration (ms)','FontName','Arial','FontSize',10);
    ylabel('M = R[S] + R[\ThetaS]','FontName','Arial','FontSize',10);
    stinrf_fig_style(ax, pretty(stype));
    if k == 1
        yl = ylim(ax);
        text(S.dur_list(1), yl(2), ' Flat mean', 'Color', c.right, ...
             'FontName','Arial','FontSize',9,'FontWeight','bold', ...
             'VerticalAlignment','top');
        text(S.dur_list(1), yl(1), ' Tb weighted', 'Color', c.ref, ...
             'FontName','Arial','FontSize',9,'FontWeight','bold', ...
             'VerticalAlignment','bottom');
    end

    ax = axes('Units','normalized','Position',[0.715 yb 0.24 0.31]); hold on;
    plot(S.dur_list, abs(S.ratio(k,:,1)), '-o', 'Color', c.right, ...
         'LineWidth',1.6,'MarkerFaceColor','w','MarkerSize',6);
    plot(S.dur_list, abs(S.ratio(k,:,2)), '-s', 'Color', c.ref, ...
         'LineWidth',1.6,'MarkerFaceColor','w','MarkerSize',6);
    fly = fly_ratio(cfg, stype);
    yline(fly,':','Color','k','LineWidth',1.1);
    text(S.dur_list(end), fly, sprintf('fly %.2f ', fly), 'FontName','Arial', ...
         'FontSize',8,'HorizontalAlignment','right','VerticalAlignment','bottom');
    set(gca,'YScale','log','XTick',S.dur_list);
    xlabel('Duration (ms)','FontName','Arial','FontSize',10);
    ylabel('|V| / |I|','FontName','Arial','FontSize',10);
    stinrf_fig_style(ax, 'Ratio');
end

exportgraphics(fig, fullfile(out_dir,'Fig_TimeReversalSummary.png'),'Resolution',600);
exportgraphics(fig, fullfile(out_dir,'Fig_TimeReversalSummary.pdf'),'ContentType','vector');
fprintf('\nSaved Fig_TimeReversalSummary in %s\n\n', out_dir);

function t = pretty(stype)
% PRETTY  Human-readable stimulus name.
    switch stype
        case 'LightEdges',       t = 'Light Edges';
        case 'UnevenSquareWave', t = 'Uneven Square Wave';
        otherwise,               t = stype;
    end
end

function r = fly_ratio(cfg, stype)
% FLY_RATIO  Ratio measured in the fly, for reference.
    switch stype
        case 'LightEdges',       r = cfg.stim.ratio_fly_LightEdges;
        case 'UnevenSquareWave', r = cfg.stim.ratio_fly_UnevenSquareWave;
        otherwise,               r = NaN;
    end
end
