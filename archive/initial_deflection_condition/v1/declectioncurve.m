% --- 初期設定 ---
clear; % ワークスペースの変数をクリア
clc;   % コマンドウィンドウをクリア
close all; % すべてのフィギュアウィンドウを閉じる

% --- パラメータ設定 ---
% 💡 ここで梁の「実際の長さ」と「最大たわみ量」を設定します
L = 18;          % [mm] 梁の実際の長さ
delta_0 = 0.05;  % [mm] グラフの最大値（最大たわみ量）

% --- 定数の定義 ---
beta_l = 4.73004;
K = 0.9825;

% --- 計算の実行 ---
% 1. 正規化された座標 (x/L) を作成 (0から1まで)
x_norm = linspace(0, 1, 500);

% 2. 実際のx軸の座標を計算 ✅
%    正規化された座標に、実際の長さLを掛ける
x_actual = x_norm * L;

% 3. 式中のβxを計算
beta_x = beta_l * x_norm;

% 4. たわみ曲線の基本形状y(x)を計算
y = cosh(beta_x) - cos(beta_x) - K * (sinh(beta_x) - sin(beta_x));

% 5. 正規化 (Normalization)
%    計算結果の形状(最大値が1のカーブ)を取り出す
y_norm = y / max(y);

% 6. 実際のたわみ曲線 w(x) を計算
%    形状(y_norm)に振幅(delta_0)を掛ける
w = delta_0 * y_norm;

% --- グラフの描画 ---
figure; % 新しいウィンドウでグラフを表示
plot(x_actual, w, 'LineWidth', 2); % ✅ x軸に実際の座標(x_actual)を使用

% グラフの体裁を整える
grid on; % グリッド線を表示
xlabel('Position x [mm]', 'FontSize', 12); % ✅ 横軸ラベルを修正
ylabel('Deflection w(x) [mm]', 'FontSize', 12);
title(['Deflection Curve (L = ' num2str(L) 'mm, Max Deflection = ' num2str(delta_0) 'mm)'], 'FontSize', 14); % ✅ タイトルを修正

% 軸の範囲を調整
xlim([0, L]);
ylim([0, delta_0 * 1.1]); % y軸の上限を最大値の1.1倍に設定