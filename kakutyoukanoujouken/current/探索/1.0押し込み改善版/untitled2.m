% 論文用の図を生成するコード
figure('Position', [50 50 1000 400], 'Color', 'w');

subplot(1,2,1);
scatter(xp_all, zeq_all, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Push-in amount x_{push} [mm]', 'FontSize', 12);
ylabel('Initial deflection z_{eq} [mm]', 'FontSize', 12);
title('Validation of theoretical formula', 'FontSize', 14);
hold on;
xp_theory = linspace(min(xp_all), max(xp_all), 100);
zeq_theory = sqrt(2*xp_theory/0.2794);
plot(xp_theory, zeq_theory, 'r-', 'LineWidth', 2);
legend('Simulation', 'Theory: z_{eq}=\sqrt{2x_p/0.2794}', 'Location', 'northwest');
grid on;
text(0.7, 1.0, sprintf('r^2 = %.4f', corr(xp_all', zeq_all')^2), 'FontSize', 12);

subplot(1,2,2);
histogram(C_eff_calc, 50, 'FaceColor', 'b', 'EdgeColor', 'k');
xlabel('Effective shape factor C_{eff}', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Distribution of C_{eff}', 'FontSize', 14);
xline(0.2794, 'r--', 'LineWidth', 2, 'Label', 'Theoretical value');
xline(mean(C_eff_calc), 'g--', 'LineWidth', 2, 'Label', sprintf('Mean = %.4f', mean(C_eff_calc)));
grid on;

print('paper_figure_validation.png', '-dpng', '-r300');
