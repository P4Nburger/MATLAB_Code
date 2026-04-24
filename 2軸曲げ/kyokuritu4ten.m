% =========================================================
% 最初の4点の(x, y)座標から根本(x=0)の曲率を計算するプログラム
% =========================================================

clear; clc; close all;

% --- 1. データの準備 ---
% x座標とy座標（画像から推測したような最初の4点のダミーデータ）
% ※ご自身の実際のデータに書き換えてください。
x = [0, 0.5, 1.0, 1.5];           % 縦方向(mm)
y = [0, -0.005, -0.019, -0.040];  % 変位(mm)

% --- 2. 曲率の計算 ---

% 【方法1】有限差分法（前進差分・2次精度）
% 条件: xが等間隔(h)であること
h = x(2) - x(1);
% 公式: y''(0) = (2*y(1) - 5*y(2) + 4*y(3) - y(4)) / h^2
curvature_diff = (2*y(1) - 5*y(2) + 4*y(3) - y(4)) / (h^2);

% 【方法2】多項式近似（カーブフィッティング）
% 4点あるので、3次関数 (y = ax^3 + bx^2 + cx + d) にフィッティング
% p = [a, b, c, d] が出力される
p = polyfit(x, y, 3);
% 2階微分 y'' = 6ax + 2b
% 根本 x=0 における2階微分（曲率）は 2b (つまり係数配列 p の2番目の要素の2倍)
curvature_poly = 2 * p(2);

% --- 3. 結果のコマンドウィンドウ表示 ---
fprintf('--- 根本(x=0)の曲率 y''''(0) の計算結果 ---\n');
fprintf('有限差分法による曲率 : %f\n', curvature_diff);
fprintf('多項式近似による曲率 : %f\n', curvature_poly);
fprintf('(※xが等間隔の場合、上記2つは数学的に一致します)\n\n');

% --- 4. プロット（視覚的な確認） ---
figure('Name', '曲率計算の確認', 'Color', 'w');
hold on; grid on;

% 4つのデータ点をプロット
plot(x, y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', '抽出した4点');

% 近似曲線のプロット（滑らかに描画するため細かいxを用意）
x_fine = linspace(min(x), max(x), 100);
y_fit = polyval(p, x_fine);
plot(x_fine, y_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', '3次多項式近似曲線');

% 根本の接線のプロット（傾き確認用: y'(0) = c）
slope = p(3); % c (xの係数)
y_tangent = slope * x_fine + y(1);
plot(x_fine, y_tangent, 'b--', 'LineWidth', 1.2, 'DisplayName', '根本の接線');

% グラフの装飾
title('根本付近のデータ点と近似曲線の確認', 'FontSize', 14);
xlabel('Longitudinal direction x [mm]', 'FontSize', 12);
ylabel('Displacement y [mm]', 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 11);

% たわみのグラフは下向きに伸びることが多いので、直感に合わせる場合Y軸を反転させることもあります
% set(gca, 'YDir', 'reverse');