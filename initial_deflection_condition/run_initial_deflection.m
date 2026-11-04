%RUN_INITIAL_DEFLECTION  フレームの押し込み量から膜の初期たわみを求める
%
%   FEA の結果で補正した理論式（式3.5）
%       x_push = (1/2) * C_eff * δ0^2      [mm]
%   を δ0 について解き、押し込み量 x_push に対する初期たわみ δ0 を計算する。
%   単位はすべて mm。
%
%   旧ファイル名: kotei_LikoruS_danmen.m（旧 kotei_LikoruS_danmen2.m もこのファイルに統合）
%   論文: 第3章（式3.5）

clear; clc; close all;

%% 1. 設定
model = 'diamond';      % 'diamond'（ひし形モデル）または 'strip'（短冊モデル）

% 押し込み量 [mm]。ベクトルで与えると一度にまとめて計算する
x_push = 0.1034;
% x_push = 0.01:0.01:0.1;   % 旧 kotei_LikoruS_danmen2.m の設定

%% 2. 有効形状定数（FEA の結果から逆算した値）
switch model
    case 'diamond'
        C_eff = 0.2795;
        % C_eff = 0.2781;   % 旧 kotei_LikoruS_danmen2.m で使っていた値
        model_name = 'Diamond Shape Model';
    case 'strip'
        C_eff = 0.2548;
        model_name = 'Strip Model';
    otherwise
        error('model は ''diamond'' か ''strip'' を指定してください。');
end

%% 3. 計算（式3.5）
delta_0 = sqrt(2 * x_push / C_eff);   % [mm] 初期たわみ

%% 4. 結果の表示
fprintf('=======================================================\n');
fprintf('<< 半経験的モデルによる予測結果 (%s) >>\n', model_name);
fprintf('使用した有効形状定数 C_eff = %.4f\n\n', C_eff);
fprintf('押し込み量[mm] | 予測たわみ量[mm]\n');
fprintf('----------------|------------------\n');
for i = 1:numel(x_push)
    fprintf('     %.4f      |      %.4f\n', x_push(i), delta_0(i));
end
fprintf('=======================================================\n');

%% 5. プロット
figure;
plot(x_push, delta_0, '-o', 'LineWidth', 2);
grid on;
xlabel('押し込み量 (x_{push}) [mm]');
ylabel('予測される最大たわみ量 (\delta_0) [mm]');
title(['Semi-Empirical Model: Push-in vs. Deflection (' model_name ')']);
xlim([0, 1.1]);
ylim([0, max(delta_0) * 1.1]);
