% =========================================================================
% 2自由度 拡張シミュレーション (速度依存型剛性モデル)
%
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
tspan = [0 1.0];
t_sound_on = 0.2;

z_target_initial = 1.0e-3; % [m] 初期たわみ

f_drive = 120;
P0_drive = 10;

% ★★★ 剛性変化パラメータ (速度依存用に再調整) ★★★
% 速度は変位より値が大きくなる(v = 2*pi*f * x)ため、係数を小さく調整
stiffness_factor = 20;       % 剛性増加係数
amplitude_sensitivity = 1.0; % 感度
amplitude_decay = 10.0;      % 減衰 (キレを良くするため少し大きく)

params = struct( ...
    'L', 30e-3, 'w1', 10e-3, 't1', 0.02e-3, 'E1', 197e9, 'rho1', 7930, ...
    'S2', 126.87e-6, 'rho2', 1240, 'k2_base', 177.0, ...
    'delta', 0.05, 'P0', P0_drive, 'freq', f_drive, ...
    'stiffness_factor', stiffness_factor, ...
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

% --- 平衡点の計算 ---
z = z_target_initial;
term = (params.k1 * z * params.L) / (4 * z * params.K_couple);
y_eq = term + (2 * z^2 / params.L);
term_couple = params.K_couple * (y_eq - 2 * z^2 / params.L);
y_natural = y_eq + term_couple / params.k2_base;
params.y_natural = y_natural;

% -------------------------------------------------------------------------
% 3. シミュレーション実行
% -------------------------------------------------------------------------
x0 = [z_target_initial; y_eq; 0; 0; 0]; 
options = odeset('RelTol',1e-6, 'AbsTol',1e-9);
[t, x] = ode45(@(t,x) equations_2DOF_VelocityBased(t, x, params, t_sound_on), tspan, x0, options);

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
hold on; xline(t_sound_on, 'r--'); hold off;
ylabel('Z 変位 [mm]'); title('膜のたわみ (Z)'); grid on;

% Y方向 (フレーム)
subplot(3,1,2);
plot(t, y*1000, 'r-', 'LineWidth', 1.5);
hold on;
xline(t_sound_on, 'r--', 'Label', 'Sound ON');
y_base = mean(y(t<t_sound_on))*1000;
yline(y_base, 'g:', 'LineWidth', 1.5, 'Label', '初期位置');
hold off;
ylabel('Y 変位 [mm]'); title('フレームの変位 (Y): 下に行けば拡張'); grid on;

% 内部変数 (剛性変化)
subplot(3,1,3);
k1_dynamic = params.k1 .* (1 + params.stiffness_factor * A);
plot(t, k1_dynamic, 'k-');
hold on; xline(t_sound_on, 'r--'); hold off;
ylabel('膜の剛性 k1 [N/m]'); 
title('剛性の変化 (速度依存)'); 
grid on; xlabel('時間 [s]');

% =========================================================================
% 2自由度 運動方程式 (速度依存型剛性変化)
% =========================================================================
function dx = equations_2DOF_VelocityBased(t, x, p, t_sound_on)
    z  = x(1);
    y  = x(2);
    dz = x(3);
    dy = x(4);
    A  = x(5);
    
    % --- 1. 剛性の動的変化 ---
    k1_curr = p.k1 * (1 + p.stiffness_factor * A);

    % --- 2. 復元力の計算 ---
    coupling_term = y - (2 * z^2 / p.L);
    
    % Z方向:
    F_restore_z = k1_curr * z - p.K_couple * coupling_term * (4 * z / p.L);
    
    % Y方向:
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
    
    % --- 5. 振幅Aの更新 (速度 dz に依存) ---
    % 静止しているとき(dz=0)は A->0 となり、剛性は元に戻る
    dA = p.amp_sens * abs(dz) - p.amp_decay * A;
    
    dx = [dz; dy; ddz; ddy; dA];
end