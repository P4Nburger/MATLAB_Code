% --- 初期設定 ---
clear; 
clc;   
close all; 

% --- パラメータ設定 ---
% 💡 ここで梁の「本来の長さ」と「目標のたわみ量」を設定します
L = 18;          % [mm] 梁の本来の長さ（まっすぐな状態）
delta_0 = 0.6404;   % [mm] 目標とする最大たわみ量

% --- 定数の定義 ---
beta_l = 4.73004;
K = 0.9825;

% --- 計算の実行 ---

% 1. たわみの「基本形状」を計算
x_norm = linspace(0, 1, 1001); 
beta_x = beta_l * x_norm;
y = cosh(beta_x) - cos(beta_x) - K * (sinh(beta_x) - sin(beta_x));
y_norm = y / max(y);

% 2. 形状定数Cを計算
dw_ds = diff(y_norm) ./ diff(x_norm);
C = trapz(x_norm(1:end-1), dw_ds.^2);

% 3. 最大たわみ量(delta_0)から押し込み量(x_push)を計算 ✅
x_push = 0.5 * C * delta_0^2;

% 4. たわみ後のパラメータを計算
L_new = L - x_push;            
x_actual = x_norm * L_new;     
w = delta_0 * y_norm;          

% --- グラフの描画 ---
figure;
plot(x_actual, w, 'LineWidth', 2);
hold on; 
plot([0, L_new], [0, 0], 'k--', 'LineWidth', 1); 
hold off;

% --- 表示の調整 ---
grid on;
axis equal; 
xlabel(['Position x [mm] (Ends at L_{new} = ' num2str(L_new, '%.3f') 'mm)'], 'FontSize', 12);
ylabel('Deflection w(x) [mm]', 'FontSize', 12);
title(['Buckled Shape (Original Length L = ' num2str(L) 'mm, Target Deflection = ' num2str(delta_0) 'mm)'], 'FontSize', 14);
legend(['Required Push-in x_{push} = ' num2str(x_push, '%.4f') ' mm'], 'Location', 'north');

% --- 結果をコマンドウィンドウに表示 ---
fprintf('===================================================\n');
fprintf('目標の最大たわみ量 δ₀ = %.4f mm を得るには、\n', delta_0);
fprintf('必要な押し込み量 x_push = %.4f mm です。\n', x_push);
fprintf('===================================================\n');