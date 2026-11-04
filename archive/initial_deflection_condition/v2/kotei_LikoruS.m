% --- 初期設定 ---
clear; 
clc;   
close all; 

% --- パラメータ設定 ---
L = 18;          % [mm] 梁の本来の長さ
x_push =1;      % [mm] フレームによる押し込み量

% --- 定数の定義 (両端固定) ---
beta_l = 4.73004;
K = 0.9825;

% --- 計算の準備 ---
L_new = L - x_push;
x_norm = linspace(0, 1, 1001);
x_actual = x_norm * L_new;

% たわみの基本形状を計算
beta_x = beta_l * x_norm;
y = cosh(beta_x) - cos(beta_x) - K * (sinh(beta_x) - sin(beta_x));
y_norm = y / max(y);

% --- ここからが新しいアプローチ --- ✅

% 1. 「誤差関数」を定義
%    入力: delta_0_guess (δ₀の推測値)
%    出力: S - L (弧長と元の長さの差)
error_func = @(delta_0_guess) ...
    (trapz(x_actual(1:end-1), sqrt(1 + (diff(delta_0_guess * y_norm) ./ diff(x_actual)).^2))) - L;

% 2. 解を探索
%    fsolveに誤差関数と、解の初期推測値を与える
initial_guess = 0.64; % 近似式で求めた値が良い出発点になる
options = optimoptions('fsolve', 'Display', 'off'); % 計算過程の表示をオフに
delta_0_exact = fsolve(error_func, initial_guess, options); % 解の探索を実行

% --- 計算結果 ---
w = delta_0_exact * y_norm;

% --- グラフの描画 ---
figure;
plot(x_actual, w, 'LineWidth', 2, 'DisplayName', ['Max Deflection \delta_0 = ' num2str(delta_0_exact, '%.4f') ' mm']);
hold on; 
plot([0, L_new], [0, 0], 'k--', 'HandleVisibility', 'off'); 
hold off;
axis equal;
grid on;
legend;
xlabel(['Position x [mm] (Ends at L_{new} = ' num2str(L_new, '%.3f') 'mm)'], 'FontSize', 12);
ylabel('Deflection w(x) [mm]', 'FontSize', 12);
title(['Buckled Shape (Exact Solution, L = ' num2str(L) 'mm, Push-in = ' num2str(x_push) 'mm)'], 'FontSize', 14);

% --- 結果の表示 ---
final_S = trapz(x_actual(1:end-1), sqrt(1 + (diff(w) ./ diff(x_actual)).^2));
fprintf('===================================================\n');
fprintf('<< 近似なしの解探索 結果 >>\n\n');
fprintf('近似式を使った場合のたわみ量 : %.4f mm\n', sqrt(2 * x_push / 4.8774));
fprintf('近似なしで求めた正確なたわみ量 : %.4f mm\n\n', delta_0_exact);
fprintf('元の長さ (L)                 : %.4f mm\n', L);
fprintf('最終的な弧長 (S)             : %.4f mm\n', final_S);
fprintf('===================================================\n');