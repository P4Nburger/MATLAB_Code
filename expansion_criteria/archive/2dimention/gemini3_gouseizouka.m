% =========================================================================
% 2自由度 拡張シミュレーション (しきい値付き・物理準拠版)
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
tspan = [0 1.0];
t_sound_on = 0.2;   % [s] 音波開始

z_target_initial = 1.0e-3; % [m] 初期たわみ

f_drive = 120;
P0_drive = 25;      % [Pa] 音圧 (0またぎを起こすため強めに)

% --- ★ 剛性変化パラメータ (しきい値付き) ★ ---
stiffness_factor = 3000;     % 剛性増加係数
stiffness_threshold = 0.02;  % ★ しきい値 (この速度を超えると剛性が増える)
amplitude_sensitivity = 3.0; 
amplitude_decay = 5.0;       

params = struct( ...
    'L', 30e-3, 'w1', 10e-3, 't1', 0.02e-3, 'E1', 197e9, 'rho1', 7930, ...
    'S2', 126.87e-6, 'rho2', 1240, 'k2_base', 177.0, ...
    'delta', 0.05, 'P0', P0_drive, 'freq', f_drive, ...
    'stiffness_factor', stiffness_factor, ...
    'threshold', stiffness_threshold, ... % パラメータに追加
    'amp_sens', amplitude_sensitivity, 'amp_decay', amplitude_decay ...
);

% -------------------------------------------------------------------------
% 2. 物理定数と平衡点の計算
% -------------------------------------------------------------------------
params.m1 = params.rho1 * params.L * params.w1 * params.t1;
params.m2 = params.rho2 * params.S2 * params.L;
params.I1 = params.w1 * params.t1^3 / 12;
params.k1 = 384 * params.E1 * params.I1 / (params.L^3);
params.K_couple = params.E1 * (params.w1 * params.t1) / params.L;
params.c1 = 2 * params.delta * sqrt(params.m1 * params.k1);
params.c2 = 2 * params.delta * sqrt(params.m2 * params.k2_base);

% --- 平衡点の計算 (ここが重要) ---
% 初期状態(A=0)での剛性 k1 を基準に平衡点を計算する
z = z_target_initial;
term = (params.k1 * z * params.L) / (4 * z * params.K_couple);
y_eq = term + (2 * z^2 / params.L);
term_couple = params.K_couple * (y_eq - 2 * z^2 / params.L);
y_natural = y_eq + term_couple / params.k2_base;
params.y_natural = y_natural;
params.z_eq = z;
params.y_eq = y_eq;

fprintf('平衡点計算完了: z=%.3f mm, y=%.3f mm\n', z*1000, y_eq*1000);

% -------------------------------------------------------------------------
% 3. シミュレーション実行
% -------------------------------------------------------------------------
% 初期条件: 平衡点からスタート (速度とAはゼロ)
x0 = [z_target_initial; y_eq; 0; 0; 0]; 

options = odeset('RelTol',1e-6, 'AbsTol',1e-9);
[t, x] = ode45(@(t,x) equations_2DOF_Threshold(t, x, params, t_sound_on), tspan, x0, options);

z = x(:,1);
y = x(:,2);
A = x(:,5);

% -------------------------------------------------------------------------
% 4. プロット
% -------------------------------------------------------------------------
figure('Position', [100, 100, 800, 900], 'Color', 'w');

% Z方向 (膜)
subplot(3,1,1);
plot(t, z*1000, 'b-');
hold on; 
xline(t_sound_on, 'r--', 'Label', 'Sound ON');
yline(0, 'k-', 'LineWidth', 1); % z=0ライン
hold off;
ylabel('Z 変位 [mm]'); title('膜のたわみ (Z): しきい値モデル'); grid on;

% Y方向 (フレーム)
subplot(3,1,2);
y_extension = (params.y_eq - y) * 1000;
plot(t, y_extension, 'r-', 'LineWidth', 1.5);
hold on;
xline(t_sound_on, 'r--');
yline(params.y_eq*1000, 'g:', 'LineWidth', 1.5, 'Label', '初期位置');
hold off;
ylabel('Y 変位 [mm]'); title('フレームの変位 (Y)'); grid on;

% 内部変数 (剛性変化)
subplot(3,1,3);
% グラフ用剛性再計算
A_effective = max(0, A - params.threshold);
k1_dynamic = params.k1 .* (1 + params.stiffness_factor * A_effective);

plot(t, k1_dynamic, 'k-', 'LineWidth', 1.5);
hold on; xline(t_sound_on, 'r--'); hold off;
ylabel('膜の剛性 k1 [N/m]'); 
title('剛性の変化 (しきい値を超えると増加)'); 
grid on; xlabel('時間 [s]');

% =========================================================================
% 2自由度 運動方程式 (しきい値付き)
% =========================================================================
function dx = equations_2DOF_Threshold(t, x, p, t_sound_on)
    z  = x(1);
    y  = x(2);
    dz = x(3);
    dy = x(4);
    A  = x(5); % 振動の激しさ
    
    % --- 1. 剛性の動的変化 (しきい値ロジック) ---
    % A が threshold を超えた分だけ剛性を増やす
    A_effective = max(0, A - p.threshold);
    k1_curr = p.k1 * (1 + p.stiffness_factor * A_effective);

    % --- 2. 復元力の計算 ---
    coupling_term = y - (2 * z^2 / p.L);
    
    % Z方向
    F_restore_z = k1_curr * z - p.K_couple * coupling_term * (4 * z / p.L);
    
    % Y方向
    F_restore_y = p.k2_base * (y - p.y_natural) + p.K_couple * coupling_term;
    
    % --- 3. 外力 ---
    if t >= t_sound_on
        F_sound = p.P0 * (p.L * p.w1) * cos(2 * pi * p.freq * t);
    else
        F_sound = 0;
    end
    
    % --- 4. 運動方程式 ---
    ddz = (F_sound - p.c1 * dz - F_restore_z) / p.m1;
    ddy = (0 - p.c2 * dy - F_restore_y) / p.m2;
    
    % --- 5. 振幅Aの更新 (速度ベース) ---
    dA = p.amp_sens * abs(dz) - p.amp_decay * A;
    
    dx = [dz; dy; ddz; ddy; dA];
end