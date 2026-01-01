% =========================================================================
% 2自由度 (2-DOF) 拡張メタマテリアル シミュレーション
% 物理モデル: ラグランジュ方程式による連成振動モデル
% 【音波ありバージョン】: 拡張現象の再現
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
% 時間設定
tspan = [0 1.0];     % [s] シミュレーション時間
t_sound_on = 0.2;    % [s] 音波をかけ始める時間

% --- 目標とする状態の設定 ---
z_target_initial = 1.0e-3; % [m] 初期たわみ量 (1mm)

% --- 音波の設定 ---
f_drive = 120;       % [Hz] 駆動周波数 (共振しそうな値)
P0_drive = 15;       % [Pa] 音圧振幅 (少し強めに設定)

% --- 物理定数 ---
params = struct();
params.L = 30e-3;      % [m] 膜の長さ
params.w1 = 10e-3;     % [m] 膜の幅
params.t1 = 0.02e-3;   % [m] 膜の厚さ
params.E1 = 197e9;     % [Pa] 膜のヤング率
params.rho1 = 7930;    % [kg/m^3] 膜の密度

params.S2 = 126.87e-6; % [m^2] フレーム断面積
params.rho2 = 1240;    % [kg/m^3] フレーム密度
params.k2 = 177.0;     % [N/m] フレームのばね定数

params.delta = 0.05;   % 減衰比
params.P0 = P0_drive;
params.freq = f_drive;

% -------------------------------------------------------------------------
% 2. 派生パラメータと平衡点の計算
% -------------------------------------------------------------------------
% 質量と剛性
params.m1 = params.rho1 * params.L * params.w1 * params.t1;
params.m2 = params.rho2 * params.S2 * params.L;
params.I1 = params.w1 * params.t1^3 / 12;
params.k1 = 384 * params.E1 * params.I1 / (params.L^3);

% ★ 連成剛性 (K_couple) ★
params.K_couple = params.E1 * (params.w1 * params.t1) / params.L;

% 減衰係数
params.c1 = 2 * params.delta * sqrt(params.m1 * params.k1);
params.c2 = 2 * params.delta * sqrt(params.m2 * params.k2);

% --- 平衡点の逆算 ---
% 1. Z方向の釣り合いから y_eq を計算
z = z_target_initial;
term = (params.k1 * z * params.L) / (4 * z * params.K_couple);
y_eq = term + (2 * z^2 / params.L);

% 2. Y方向の釣り合いから y_natural を計算
term_couple = params.K_couple * (y_eq - 2 * z^2 / params.L);
y_natural = y_eq + term_couple / params.k2;

params.y_natural = y_natural;

fprintf('--- 設定完了 ---\n');
fprintf('初期たわみ (Z): %.3f mm\n', z*1000);
fprintf('初期収縮量 (Y): %.3f mm\n', y_eq*1000);
fprintf('音圧 (P0): %.1f Pa\n', P0_drive);

% -------------------------------------------------------------------------
% 3. シミュレーション実行
% -------------------------------------------------------------------------
x0 = [z_target_initial; y_eq; 0; 0]; 

fprintf('シミュレーション実行中...\n');
options = odeset('RelTol',1e-6, 'AbsTol',1e-9);
[t, x] = ode45(@(t,x) equations_2DOF(t, x, params, t_sound_on), tspan, x0, options);

z_res = x(:,1);
y_res = x(:,2);

% -------------------------------------------------------------------------
% 4. プロット
% -------------------------------------------------------------------------
figure('Position', [100, 100, 800, 800], 'Color', 'w');

% (1) Z方向 (膜)
subplot(2,1,1);
plot(t, z_res*1000, 'b-', 'LineWidth', 1.0);
hold on;
xline(t_sound_on, 'r--', 'LineWidth', 1.5, 'Label', 'Sound ON');
yline(z_target_initial*1000, 'g:', 'LineWidth', 1.5, 'Label', '初期位置');
hold off;
ylabel('Z方向変位 [mm]', 'FontSize', 12);
title('膜のたわみ (Z)', 'FontSize', 14);
grid on;

% (2) Y方向 (フレーム)
subplot(2,1,2);
plot(t, y_res*1000, 'r-', 'LineWidth', 1.0);
hold on;
xline(t_sound_on, 'r--', 'LineWidth', 1.5, 'Label', 'Sound ON');
y_base = y_eq * 1000;
yline(y_base, 'g:', 'LineWidth', 1.5, 'Label', '初期位置');
hold off;
ylabel('Y方向変位 (収縮量) [mm]', 'FontSize', 12);
title('フレームの変位 (Y): 下に行くほど拡張', 'FontSize', 14);
grid on;
xlabel('時間 [s]', 'FontSize', 12);

% 結果の判定
y_after = mean(y_res(t > 0.5));
fprintf('---------------------------------------------------\n');
fprintf('初期Y位置: %.4f mm\n', y_eq*1000);
fprintf('振動中Y平均: %.4f mm\n', y_after*1000);
if y_after < y_eq
    fprintf('判定: 成功 (拡張しました)\n');
else
    fprintf('判定: 収縮または変化なし\n');
end
fprintf('---------------------------------------------------\n');


% =========================================================================
% 2自由度 運動方程式 (音波あり)
% =========================================================================
function dx = equations_2DOF(t, x, p, t_sound_on)
    % 状態変数の展開
    z  = x(1);
    y  = x(2);
    dz = x(3);
    dy = x(4);
    
    % --- 力の計算 ---
    
    % 1. 膜の幾何学的連成項 (膜の伸び量)
    coupling_term = y - (2 * z^2 / p.L);
    
    % 2. 復元力
    % Z方向
    F_restore_z = p.k1 * z - p.K_couple * coupling_term * (4 * z / p.L);
    
    % Y方向
    F_restore_y = p.k2 * (y - p.y_natural) + p.K_couple * coupling_term;
    
    % 3. 外力 (音波)
    if t >= t_sound_on
        % 面積 p.L * p.w1 に圧力 p.P0 がかかる
        F_sound = p.P0 * (p.L * p.w1) * cos(2 * pi * p.freq * t);
    else
        F_sound = 0;
    end
    
    % 4. 減衰力
    F_damping_z = p.c1 * dz;
    F_damping_y = p.c2 * dy;
    
    % --- 運動方程式 (F = ma) ---
    ddz = (F_sound - F_damping_z - F_restore_z) / p.m1;
    ddy = (0 - F_damping_y - F_restore_y) / p.m2;
    
    dx = [dz; dy; ddz; ddy];
end