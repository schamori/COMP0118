clc
clear all
close all

T2 = [5, 10, 15, 20];
TE = [1:1.375:16.5]';
s0 = [155, 255, 355, 455];
Nrow = 32;
Ncol = 32;

quadrant_labels = {'Q1: a=155, T2*=5', 'Q2: a=255, T2*=10', ...
                   'Q3: a=355, T2*=15', 'Q4: a=455, T2*=20'};
sigma2_vals  = [0, 4];
noise_labels = {'No Noise (\sigma^2=0)', 'Gaussian Noise (\sigma^2=4)'};
noise_short  = {'No Noise', 'Noisy (\sigma^2=4)'};

for uu = 1:length(T2)
    Phantom_WO_NoiseTemp{uu} = createPhantoms('exp', TE, T2(uu), s0(uu), Nrow, Ncol);
end
Phantom_WO_Noise = [Phantom_WO_NoiseTemp{1}, Phantom_WO_NoiseTemp{2}; ...
                    Phantom_WO_NoiseTemp{3}, Phantom_WO_NoiseTemp{4}];
[Nrow_, Ncol_, bands] = size(Phantom_WO_Noise);

S0_Image = [s0(1)*ones(Nrow,Ncol), s0(2)*ones(Nrow,Ncol); ...
            s0(3)*ones(Nrow,Ncol), s0(4)*ones(Nrow,Ncol)];
T2_Image = [T2(1)*ones(Nrow,Ncol), T2(2)*ones(Nrow,Ncol); ...
            T2(3)*ones(Nrow,Ncol), T2(4)*ones(Nrow,Ncol)];

Y_all = cell(2,1);
for nn = 1:2
    Y_all{nn} = Phantom_WO_Noise + sqrt(sigma2_vals(nn)) * randn(Nrow_, Ncol_, bands);
end

colors_q = lines(4);

figure(1)
set(gcf,'Units','normalized','Position',[0.02 0.5 0.96 0.45])
sgtitle('Task 4/5 — TE vs Mean Signal Intensity per Quadrant', ...
        'FontSize', 13, 'FontWeight', 'bold')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

for nn = 1:2
    Y = Y_all{nn};
    nexttile
    hold on
    for q = 1:4
        row_idx = (q > 2)*Nrow + (1:Nrow);
        col_idx = (mod(q-1,2))*Ncol + (1:Ncol);
        q_pixels = reshape(Y(row_idx, col_idx, :), [], bands);
        q_mean   = mean(q_pixels, 1);
        plot(TE, q_mean, '-o', 'Color', colors_q(q,:), 'LineWidth', 2, ...
             'MarkerSize', 5, 'DisplayName', quadrant_labels{q})
    end
    hold off
    xlabel('TE (ms)', 'FontSize', 11)
    ylabel('Mean Signal Intensity', 'FontSize', 11)
    title(noise_labels{nn}, 'FontSize', 12)
    legend('Location', 'northeast', 'FontSize', 9)
    grid on; box on
end

lambdaA = 1;
lambdaR = 1;

a_maps   = cell(2,1);
T2_maps  = cell(2,1);
a_means  = zeros(2,4);
T2_means = zeros(2,4);

for nn = 1:2
    fprintf('\n========== %s ==========\n', noise_labels{nn});
    yReshaped = reshape(Y_all{nn}, Nrow_*Ncol_, bands)';
    [a, r, ~, ~] = relaxationEst(yReshaped, TE, Nrow_, Ncol_, lambdaA, lambdaR);

    a_map  = reshape(a,    Nrow_, Ncol_);
    T2_map = reshape(1./r, Nrow_, Ncol_);
    a_maps{nn}  = a_map;
    T2_maps{nn} = T2_map;

    a_means(nn,:)  = [mean(a_map(1:32,   1:32),  'all'), mean(a_map(1:32,   33:end),'all'), ...
                      mean(a_map(33:end, 1:32),  'all'), mean(a_map(33:end, 33:end),'all')];
    T2_means(nn,:) = [mean(T2_map(1:32,  1:32),  'all'), mean(T2_map(1:32,  33:end),'all'), ...
                      mean(T2_map(33:end,1:32),  'all'), mean(T2_map(33:end,33:end),'all')];

    fprintf('  Quadrant  |  GT a   Est a   Err%%  |  GT T2*  Est T2*  Err%%\n');
    fprintf('  ----------+----------------------+-------------------------\n');
    for q = 1:4
        fprintf('  Q%d        |  %5.1f  %6.2f  %5.1f%% |  %5.1f    %6.2f   %5.1f%%\n', q, ...
            s0(q), a_means(nn,q),  100*(a_means(nn,q) - s0(q))/s0(q), ...
            T2(q),  T2_means(nn,q), 100*(T2_means(nn,q) - T2(q))/T2(q))
    end
end

figure(2)
set(gcf,'Units','normalized','Position',[0.02 0.5 0.96 0.45])
sgtitle('Task 5 — Estimated a maps', 'FontSize', 13, 'FontWeight', 'bold')
tiledlayout(1,4,'TileSpacing','compact','Padding','compact')
clim_a = [0, max(s0)*1.05];

nexttile; imagesc(S0_Image); axis image off; clim(clim_a); colorbar
title('Ground Truth a_0', 'FontSize', 11)

nexttile; imagesc(a_maps{1}); axis image off; clim(clim_a); colorbar
title('Estimated a  (no noise)', 'FontSize', 11)

nexttile; imagesc(a_maps{2}); axis image off; clim(clim_a); colorbar
title('Estimated a  (\sigma^2=4)', 'FontSize', 11)

nexttile; imagesc(abs(a_maps{2} - a_maps{1})); axis image off; colorbar; colormap(gca, hot)
title('|Noisy - Clean|  difference', 'FontSize', 11)

figure(3)
set(gcf,'Units','normalized','Position',[0.02 0.0 0.96 0.45])
sgtitle('Task 5 — Estimated T2* maps', 'FontSize', 13, 'FontWeight', 'bold')
tiledlayout(1,4,'TileSpacing','compact','Padding','compact')
clim_T2 = [0, max(T2)*1.15];

nexttile; imagesc(T2_Image); axis image off; clim(clim_T2); colorbar; colormap(gca, parula)
title('Ground Truth T2*', 'FontSize', 11)

nexttile; imagesc(T2_maps{1}); axis image off; clim(clim_T2); colorbar; colormap(gca, parula)
title('Estimated T2*  (no noise)', 'FontSize', 11)

nexttile; imagesc(T2_maps{2}); axis image off; clim(clim_T2); colorbar; colormap(gca, parula)
title('Estimated T2*  (\sigma^2=4)', 'FontSize', 11)

nexttile; imagesc(abs(T2_maps{2} - T2_maps{1})); axis image off; colorbar; colormap(gca, hot)
title('|Noisy - Clean|  difference', 'FontSize', 11)

figure(4)
set(gcf,'Units','normalized','Position',[0.02 0.5 0.46 0.45])
bar_data_a = [s0(:), a_means(1,:)', a_means(2,:)'];
b = bar(bar_data_a, 'grouped');
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.2 0.75 0.3];
b(3).FaceColor = [0.9 0.3 0.2];
hold on
for q = 1:4
    for kk = 2:3
        x_pos = b(kk).XEndPoints(q);
        err = 100*(bar_data_a(q,kk) - s0(q)) / s0(q);
        text(x_pos, bar_data_a(q,kk) + 3, sprintf('%.1f%%', err), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.2 0.2 0.2])
    end
end
hold off
set(gca, 'XTickLabel', {'Q1','Q2','Q3','Q4'}, 'FontSize', 11)
ylabel('Mean Signal Intensity a', 'FontSize', 12)
title('Mean a — GT vs Estimated', 'FontSize', 13, 'FontWeight', 'bold')
legend({'Ground Truth','No Noise','Noisy (\sigma^2=4)'}, 'Location','northwest','FontSize',10)
ylim([0, max(s0)*1.15]); grid on; box on

figure(5)
set(gcf,'Units','normalized','Position',[0.5 0.5 0.46 0.45])
bar_data_T2 = [T2(:), T2_means(1,:)', T2_means(2,:)'];
b2 = bar(bar_data_T2, 'grouped');
b2(1).FaceColor = [0.2 0.4 0.8];
b2(2).FaceColor = [0.2 0.75 0.3];
b2(3).FaceColor = [0.9 0.3 0.2];
hold on
for q = 1:4
    for kk = 2:3
        x_pos = b2(kk).XEndPoints(q);
        err = 100*(bar_data_T2(q,kk) - T2(q)) / T2(q);
        text(x_pos, bar_data_T2(q,kk) + 0.2, sprintf('%.1f%%', err), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.2 0.2 0.2])
    end
end
hold off
set(gca, 'XTickLabel', {'Q1','Q2','Q3','Q4'}, 'FontSize', 11)
ylabel('Mean T2* (ms)', 'FontSize', 12)
title('Mean T2* — GT vs Estimated', 'FontSize', 13, 'FontWeight', 'bold')
legend({'Ground Truth','No Noise','Noisy (\sigma^2=4)'}, 'Location','northwest','FontSize',10)
ylim([0, max(T2)*1.2]); grid on; box on

figure(6)
set(gcf,'Units','normalized','Position',[0.02 0.0 0.96 0.45])
sgtitle('Task 5 — Fitted decay curves: GT vs estimated', ...
        'FontSize', 13, 'FontWeight', 'bold')
tiledlayout(1,4,'TileSpacing','compact','Padding','compact')

TE_fine  = linspace(TE(1), TE(end), 200);
clrs     = lines(3);

for q = 1:4
    nexttile
    plot(TE_fine, s0(q)*exp(-TE_fine/T2(q)), 'k-', 'LineWidth', 2.5, ...
         'DisplayName', 'Ground Truth'); hold on

    for nn = 1:2
        est = a_means(nn,q) * exp(-TE_fine / T2_means(nn,q));
        plot(TE_fine, est, '--', 'Color', clrs(nn+1,:), 'LineWidth', 2, ...
             'DisplayName', noise_short{nn})
    end

    row_idx = (q > 2)*Nrow + (1:Nrow);
    col_idx = (mod(q-1,2))*Ncol + (1:Ncol);
    q_pixels = reshape(Phantom_WO_Noise(row_idx, col_idx, :), [], bands);
    scatter(TE, mean(q_pixels,1), 30, 'k', 'filled', 'HandleVisibility','off')

    hold off
    xlabel('TE (ms)', 'FontSize', 10); ylabel('Signal Intensity', 'FontSize', 10)
    title(quadrant_labels{q}, 'FontSize', 11)
    legend('Location','northeast','FontSize',8)
    grid on; box on
end
