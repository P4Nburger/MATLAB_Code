% =========================================================================
% 4変数探索結果の確認（修正版）
% =========================================================================
clear; clc;

% データ読み込み
load('results_20260119_145358\zeq_extended_search.mat');

fprintf('========================================\n');
fprintf('   4-Variable Search Results Summary\n');
fprintf('========================================\n\n');

% --- 探索範囲の確認（データから逆算） ---
all_xp = [res.x_push_m] * 1000;
all_L = [res.L_m] * 1000;
all_ws = [res.w1_stiff_m] * 1000;
all_wm = [res.w1_mass_m] * 1000;

xp_unique = unique(all_xp);
L_unique = unique(all_L);
ws_unique = unique(all_ws);
wm_unique = unique(all_wm);

fprintf('【探索範囲】\n');
fprintf('x_push:   %.3f ~ %.3f mm (%d values)\n', ...
    min(xp_unique), max(xp_unique), length(xp_unique));
fprintf('L:        %.1f ~ %.1f mm (%d values)\n', ...
    min(L_unique), max(L_unique), length(L_unique));
fprintf('w1_stiff: %.2f ~ %.1f mm (%d values)\n', ...
    min(ws_unique), max(ws_unique), length(ws_unique));
fprintf('w1_mass:  %.1f ~ %.1f mm (%d values)\n', ...
    min(wm_unique), max(wm_unique), length(wm_unique));

n_total = length(res);
fprintf('\n総組み合わせ数: %d\n', n_total);
fprintf('理論値: %d × %d × %d × %d = %d\n', ...
    length(xp_unique), length(L_unique), length(ws_unique), length(wm_unique), ...
    length(xp_unique)*length(L_unique)*length(ws_unique)*length(wm_unique));

% --- 収束状況 ---
cv = [res.converged];
n_converged = sum(cv);
n_failed = n_total - n_converged;

fprintf('\n【収束状況】\n');
fprintf('収束成功: %d / %d (%.1f%%)\n', n_converged, n_total, 100*n_converged/n_total);
fprintf('収束失敗: %d / %d (%.1f%%)\n\n', n_failed, n_total, 100*n_failed/n_total);

% --- 収束した設計の統計 ---
r = res(cv);
xp_all = [r.x_push_m] * 1000;
zeq_all = [r.z_eq_mm];
L_all = [r.L_m] * 1000;
ws_all = [r.w1_stiff_m] * 1000;
wm_all = [r.w1_mass_m] * 1000;
f_all = [r.f_res_Hz];
Az_all = [r.A_z_mm];

fprintf('【収束した設計の統計】\n');
fprintf('z_eq:     %.3f ~ %.3f mm (mean: %.3f, std: %.3f)\n', ...
    min(zeq_all), max(zeq_all), mean(zeq_all), std(zeq_all));
fprintf('x_push:   %.3f ~ %.3f mm (mean: %.3f, std: %.3f)\n', ...
    min(xp_all), max(xp_all), mean(xp_all), std(xp_all));
fprintf('A_z:      %.3f ~ %.3f mm (mean: %.3f, std: %.3f)\n', ...
    min(Az_all), max(Az_all), mean(Az_all), std(Az_all));
fprintf('f_res:    %.1f ~ %.1f Hz (mean: %.1f, std: %.1f)\n', ...
    min(f_all), max(f_all), mean(f_all), std(f_all));
fprintf('L:        %.1f ~ %.1f mm\n', min(L_all), max(L_all));
fprintf('w_stiff:  %.2f ~ %.1f mm\n', min(ws_all), max(ws_all));
fprintf('w_mass:   %.1f ~ %.1f mm\n\n', min(wm_all), max(wm_all));

% --- C_effの確認 ---
C_eff_calc = 2 * xp_all ./ (zeq_all.^2);
fprintf('【C_eff 検証】\n');
fprintf('C_eff = 2*xp/zeq^2\n');
fprintf('Mean: %.4f\n', mean(C_eff_calc));
fprintf('Std:  %.4f (%.1f%%)\n', std(C_eff_calc), 100*std(C_eff_calc)/mean(C_eff_calc));
fprintf('Min:  %.4f\n', min(C_eff_calc));
fprintf('Max:  %.4f\n\n', max(C_eff_calc));

if std(C_eff_calc)/mean(C_eff_calc) < 0.05
    fprintf('→ C_effはほぼ一定（変動 < 5%%）\n');
    fprintf('→ z_eq ≈ sqrt(2*xp/%.4f) が成立\n\n', mean(C_eff_calc));
else
    fprintf('→ C_effは設計パラメータに依存（変動 > 5%%）\n');
    fprintf('→ さらなる分析が必要\n\n');
end

% --- Top 20 by z_eq ---
[zeq_sorted, idx] = sort(zeq_all, 'descend');
top20_idx = idx(1:min(20, length(idx)));

fprintf('【Top 20 by z_eq】\n');
fprintf('Rank | z_eq [mm] | xp [mm] | L [mm] | ws [mm] | wm [mm] | f [Hz] | A_z [mm] | A_z/z_eq\n');
fprintf('-----|-----------|---------|--------|---------|---------|--------|----------|----------\n');
for i = 1:min(20, length(top20_idx))
    d = r(top20_idx(i));
    ratio = d.A_z_mm / d.z_eq_mm;
    fprintf('%4d | %9.3f | %7.3f | %6.1f | %7.2f | %7.1f | %6.1f | %8.3f | %8.2f\n', ...
        i, d.z_eq_mm, d.x_push_m*1000, d.L_m*1000, ...
        d.w1_stiff_m*1000, d.w1_mass_m*1000, d.f_res_Hz, d.A_z_mm, ratio);
end

% --- パラメータの傾向分析 ---
fprintf('\n【Top 20の傾向】\n');
top20 = r(top20_idx);
top20_xp = [top20.x_push_m]*1000;
top20_L = [top20.L_m]*1000;
top20_ws = [top20.w1_stiff_m]*1000;
top20_wm = [top20.w1_mass_m]*1000;

fprintf('x_push:   %.3f ~ %.3f mm (平均: %.3f, 最頻値: %.3f)\n', ...
    min(top20_xp), max(top20_xp), mean(top20_xp), mode(top20_xp));
fprintf('L:        %.1f ~ %.1f mm (平均: %.1f, 最頻値: %.1f)\n', ...
    min(top20_L), max(top20_L), mean(top20_L), mode(top20_L));
fprintf('w_stiff:  %.2f ~ %.1f mm (平均: %.2f, 最頻値: %.2f)\n', ...
    min(top20_ws), max(top20_ws), mean(top20_ws), mode(top20_ws));
fprintf('w_mass:   %.1f ~ %.1f mm (平均: %.1f, 最頻値: %.1f)\n', ...
    min(top20_wm), max(top20_wm), mean(top20_wm), mode(top20_wm));

% 各パラメータの出現頻度
fprintf('\n【Top 20での出現頻度】\n');
fprintf('x_push = %.3f mm: %d回\n', mode(top20_xp), sum(top20_xp == mode(top20_xp)));
fprintf('L = %.1f mm: %d回\n', mode(top20_L), sum(top20_L == mode(top20_L)));
fprintf('w_stiff = %.2f mm: %d回\n', mode(top20_ws), sum(top20_ws == mode(top20_ws)));
fprintf('w_mass = %.1f mm: %d回\n', mode(top20_wm), sum(top20_wm == mode(top20_wm)));

fprintf('\n========================================\n');

% --- 視覚化 ---
figure('Position', [50 50 1400 900], 'Color', 'w');

% (1) z_eq分布
subplot(2,3,1);
histogram(zeq_all, 30, 'FaceColor', 'b', 'EdgeColor', 'k');
hold on;
xline(mean(zeq_all), 'r--', 'LineWidth', 2, 'Label', 'Mean');
xline(zeq_sorted(20), 'g--', 'LineWidth', 2, 'Label', 'Top 20');
hold off;
xlabel('z_{eq} [mm]'); ylabel('Count');
title('Distribution of z_{eq}');
legend('Location', 'best');
grid on;

% (2) x_push vs z_eq
subplot(2,3,2);
scatter(xp_all, zeq_all, 20, zeq_all, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
% Top 20をハイライト
scatter(xp_all(top20_idx(1:10)), zeq_all(top20_idx(1:10)), 100, 'r', 'LineWidth', 2);
hold off;
xlabel('x_{push} [mm]'); ylabel('z_{eq} [mm]');
title('x_{push} vs z_{eq}');
colorbar; grid on;

% (3) L vs z_eq
subplot(2,3,3);
scatter(L_all, zeq_all, 20, zeq_all, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(L_all(top20_idx(1:10)), zeq_all(top20_idx(1:10)), 100, 'r', 'LineWidth', 2);
hold off;
xlabel('L [mm]'); ylabel('z_{eq} [mm]');
title('L vs z_{eq}');
colorbar; grid on;

% (4) w_stiff vs z_eq
subplot(2,3,4);
scatter(ws_all, zeq_all, 20, zeq_all, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(ws_all(top20_idx(1:10)), zeq_all(top20_idx(1:10)), 100, 'r', 'LineWidth', 2);
hold off;
xlabel('w_{stiff} [mm]'); ylabel('z_{eq} [mm]');
title('w_{stiff} vs z_{eq}');
colorbar; grid on;

% (5) w_mass vs z_eq
subplot(2,3,5);
scatter(wm_all, zeq_all, 20, zeq_all, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(wm_all(top20_idx(1:10)), zeq_all(top20_idx(1:10)), 100, 'r', 'LineWidth', 2);
hold off;
xlabel('w_{mass} [mm]'); ylabel('z_{eq} [mm]');
title('w_{mass} vs z_{eq}');
colorbar; grid on;

% (6) C_eff distribution
subplot(2,3,6);
histogram(C_eff_calc, 30, 'FaceColor', 'm', 'EdgeColor', 'k');
hold on;
xline(mean(C_eff_calc), 'r--', 'LineWidth', 2, 'Label', sprintf('Mean=%.4f', mean(C_eff_calc)));
xline(0.2794, 'g--', 'LineWidth', 2, 'Label', 'Theory=0.2794');
hold off;
xlabel('C_{eff}'); ylabel('Count');
title('Distribution of C_{eff}');
legend('Location', 'best');
grid on;

sgtitle('4-Variable Search Results Overview', 'FontSize', 14, 'FontWeight', 'bold');
