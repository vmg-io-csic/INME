%% plot_flash.m
%
% Generates a 3x5 figure from the results of run_flash_2D.m.

clear; clc; close all;
load('results_flash_2D.mat');

n_keys = numel(key_frames);
dpi_val = N / dov;

% Spatial axes in degrees
x_deg = ((1:N) - (N+1)/2) / dpi_val;
y_deg = ((1:N) - (N+1)/2) / dpi_val;

% Arrow positions in degrees
arr_x_deg = (arr_x_px - cx) / dpi_val;
arr_y_deg = (arr_y_px - cy) / dpi_val;

% Disc outline (circle in degrees)
th_circ = linspace(0, 2*pi, 100);
circ_x = disc_radius_deg * cos(th_circ);
circ_y = disc_radius_deg * sin(th_circ);

% Global magnitude range for consistent colormap
mag_max = max(mag_t(:));
if mag_max == 0, mag_max = 1; end

%% FIGURE
fig = figure('Position', [50 50 1400 700], 'Color', 'w');

for k = 1:n_keys
    fi = key_frames(k);

    % --- Row 1: Stimulus ---
    ax1 = subplot(3, n_keys, k);
    stim_frame = double(TCB_sq(:,:,fi));
    imagesc(x_deg, y_deg, stim_frame);
    hold on; plot(circ_x, circ_y, 'g', 'LineWidth', 1.5); hold off;
    colormap(ax1, gray); caxis([0 200]);
    axis image; set(gca, 'YDir', 'normal');
    title(key_labels{k}, 'FontSize', 10);
    if k == 1, ylabel('Stimulus', 'FontWeight', 'bold'); end
    if k == n_keys
        cb = colorbar; cb.Label.String = 'Luminance (a.u.)';
    end

    % --- Row 2: Magnitude heatmap ---
    ax2 = subplot(3, n_keys, n_keys + k);
    fm = mag_t(:,:,fi);
    imagesc(x_deg, y_deg, fm);
    hold on; plot(circ_x, circ_y, 'g', 'LineWidth', 1.5); hold off;
    colormap(ax2, hot); caxis([0 mag_max]);
    axis image; set(gca, 'YDir', 'normal');
    if k == 1, ylabel('Magnitude', 'FontWeight', 'bold'); end
    if k == n_keys
        cb = colorbar; cb.Label.String = 'Magnitude (a.u.)';
    end
    % Ring mag number
    ring_mag_k = mean(fm(mask_ring));
    text(0.95*max(x_deg), 0.90*max(y_deg), sprintf('%.1f', ring_mag_k), ...
        'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right');

    % --- Row 3: Direction arrows ---
    ax3 = subplot(3, n_keys, 2*n_keys + k);
    plot(circ_x, circ_y, 'g-', 'LineWidth', 1.5); hold on;
    axis equal;
    lim = disc_radius_deg * 1.8;
    xlim([-lim lim]); ylim([-lim lim]);
    set(gca, 'Color', 'w');

    oh_frame = O_H_t(:,:,fi);
    ov_frame = O_V_t(:,:,fi);

    for a = 1:n_arrows
        ah = oh_frame(arr_y_px(a), arr_x_px(a));
        av = ov_frame(arr_y_px(a), arr_x_px(a));
        am = sqrt(ah^2 + av^2);

        if am > 0.01 * ring_mag_k && ring_mag_k > 0.01 * mag_max
            % Normalise arrow length
            scale = 0.15 * disc_radius_deg;
            quiver(arr_x_deg(a), arr_y_deg(a), ...
                ah/am * scale, av/am * scale, 0, ...
                'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.8);
        else
            plot(arr_x_deg(a), arr_y_deg(a), 'g.', 'MarkerSize', 8);
        end
    end
    hold off;
    if k == 1, ylabel('Direction', 'FontWeight', 'bold'); end
    set(gca, 'XTick', [], 'YTick', []);
end

% Bottom annotation
annotation('textbox', [0.1, 0.01, 0.8, 0.03], ...
    'String', sprintf('Nm = %d ms  |  Disc: %.1f deg  |  Number: mean magnitude on ring  |  %d arrows at uniform angular positions (dots = no response)', ...
    model.Nm, disc_radius_deg, n_arrows), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);

%% SAVE FIGURE
print(fig, 'fig_flash_results', '-dpng', '-r300');
fprintf('Saved: fig_flash_results.png\n');
