% --- 初期設定 ---
clear; 
clc;   
close all; 

% --- パラメータ設定 ---
% 💡 ここで梁の「本来の長さ」と「押し込み量」を設定します
L = 18;          % [mm] 梁の本来の長さ（まっすぐな状態）
x_push = 1.6;    % [mm] フレームによる押し込み量

% --- 計算の実行 ---

% 1. たわみの「基本形状」を計算 (両端回転支持の場合)
%    回転支持の形状はサイン波になります。
x_norm = linspace(0, 1, 1001); 
y_norm = sin(pi * x_norm);

% 2. 形状定数Cを計算
%    C = ∫[0→1] (w_norm'(s))^2 ds
dw_ds = diff(y_norm) ./ diff(x_norm); 
C = trapz(x_norm(1:end-1), dw_ds.^2);

% 3. 押し込み量(x_push)から最大たわみ量(delta_0)を計算 ✅
delta_0 = sqrt(2 * x_push / C);

% 4. たわみ後のパラメータを計算
L_new = L - x_push;            % たわみ後の両端間距離
x_actual = x_norm * L_new;     % 実際のx軸を作成
w = delta_0 * y_norm;          % 実際のたわみ曲線を計算

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
title(['Buckled Shape (Original Length L = ' num2str(L) 'mm, Push-in = ' num2str(x_push) 'mm)'], 'FontSize', 14);
legend(['Max Deflection \delta_0 = ' num2str(delta_0, '%.4f') ' mm'], 'Location', 'north');

% --- 結果をコマンドウィンドウに表示 ---
fprintf('===================================================\n');
fprintf('<< 両端回転支持モデル >>\n');
fprintf('押し込み量 x_push = %.4f mm の時、\n', x_push);
fprintf('最大たわみ量 δ₀ = %.4f mm となります。\n', delta_0);
fprintf('形状定数 C = %.4f\n', C);
fprintf('===================================================\n');