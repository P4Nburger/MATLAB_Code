% --- 初期設定 ---
clear; 
clc;   
close all; 

% --- パラメータ設定 ---
% 💡 ここで梁の「本来の長さ」と「押し込み量」を設定します
L = 18;          % [mm] 梁の本来の長さ（まっすぐな状態）
x_push = 1;      % [mm] フレームによる押し込み量

% --- 定数の定義 (両端固定) ---
beta_l = 4.73004;
K = 0.9825;

% --- 計算の実行 ---

% 1. たわみの「基本形状」を計算 (両端固定の場合)
x_norm = linspace(0, 1, 1001); 
beta_x = beta_l * x_norm;
y = cosh(beta_x) - cos(beta_x) - K * (sinh(beta_x) - sin(beta_x));
y_norm = y / max(y);

% 2. 形状定数Cを計算
dw_ds = diff(y_norm) ./ diff(x_norm);
C = trapz(x_norm(1:end-1), dw_ds.^2);

% 3. 押し込み量(x_push)から最大たわみ量(delta_0)を計算
delta_0 = sqrt(2 * x_push / C);

% 4. たわみ後のパラメータを計算
L_new = L - x_push;            
x_actual = x_norm * L_new;     
w = delta_0 * y_norm;          

% --- 弧長の計算 (ここからが追加部分) --- ✅

% 5. たわみ曲線の微分 dw/dx を数値的に計算
dw_dx = diff(w) ./ diff(x_actual);

% 6. 積分の中身 √(1 + (dw/dx)²) を計算
integrand = sqrt(1 + dw_dx.^2);

% 7. 数値積分して弧長 S を求める
%    diffで要素数が1つ減るため、x_actualもそれに合わせる
S = trapz(x_actual(1:end-1), integrand);


% --- グラフの描画 ---
figure;
plot(x_actual, w, 'LineWidth', 2, 'DisplayName', ['Max Deflection \delta_0 = ' num2str(delta_0, '%.4f') ' mm']);
hold on; 
plot([0, L_new], [0, 0], 'k--', 'HandleVisibility', 'off');
hold off;
grid on;
legend;
xlabel(['Position x [mm] (Ends at L_{new} = ' num2str(L_new, '%.3f') 'mm)'], 'FontSize', 12);
ylabel('Deflection w(x) [mm]', 'FontSize', 12);
title(['Buckled Shape (Original Length L = ' num2str(L) 'mm, Push-in = ' num2str(x_push) 'mm)'], 'FontSize', 14);
ylim([-6, 7]); % グラフの見た目を画像に合わせる


% --- 計算結果の表示 (ここからが追加部分) --- ✅
fprintf('===================================================\n');
fprintf('<< 弧長の検証計算結果 >>\n\n');
fprintf('入力した元の長さ (L)      : %.4f mm\n', L);
fprintf('計算されたたわみ後の弧長 (S) : %.4f mm\n\n', S);
fprintf('両者の差 (L - S)            : %.6f mm\n', L -  S);
fprintf('===================================================\n');