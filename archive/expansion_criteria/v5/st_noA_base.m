% =========================================================================
% 2自由度 シミュレーション (自立維持 + 剛性一定 + 周波数スイープ)
% ・初期たわみ: フレームの内力(y_natural)で自立維持 
% ・剛性変化: ★なし (k1は常に一定)★
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
% --- スイープ設定 ---
freq_range = 20:1:30;  % [Hz]
tspan_sweep = [0 2.0];  % [s]

% --- 音波設定 ---
P0_drive = 2;          % [Pa]

% --- 本番シミュレーション設定 ---
tspan_main = [0 5.0];   % [s]
t_sound_on = 0.2;       % [s]


% --- 初期たわみ設定 (実験値準拠) ---
x_push = 0.1e-3;        % [m] フレーム押し込み量
C_eff = 0.2794;         % 有効形状定数
z_target_initial_mm = sqrt(2 * (x_push * 1000) / C_eff);
z_target_initial = z_target_initial_mm / 1000; % [m]


params = struct( ...
    'L', 30e-3, 'w1', 10e-3, 't1', 0.02e-3, 'E1', 197e9, 'rho1', 7930, ...
    'S2', 126.87e-6, 'rho2', 1240, 'k2_base', 177.0, ...
    'delta', 0.05, 'P0', P0_drive, 'freq', 0, ... 
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

% 幾何学係数 Gamma
params.Gamma = (C_eff / 2) * 1000; 

% --- 平衡点と自然長の逆算（自立維持のための計算） ---
z = z_target_initial;
% 釣り合い状態での幾何学的収縮量
y_eq = params.Gamma * z^2;

% 「フレームが縮もうとする力」と「膜の張力」を釣り合わせるための y_natural を計算
% 1. Z方向の釣り合いから、必要な張力項の値を逆算
%    F_restore_z = k1*z - K_couple * term * (2*G*z) = 0
%    -> term = k1*z / (K_couple * 2*G*z) = k1 / (2*G*K_couple)
geom_term_needed = params.k1 / (2 * params.Gamma * params.K_couple);

% 2. その張力を発生させるための y_natural をY方向の釣り合いから逆算
%    F_restore_y = k2*(y - y_nat) + K_couple * term = 0
%    -> k2*(y - y_nat) = -K_couple * term
%    -> y - y_nat = - (K_couple * term) / k2
%    -> y_nat = y + (K_couple * term) / k2
y_natural = y_eq + (params.K_couple * geom_term_needed) / params.k2_base;

params.y_eq = y_eq;
params.z_eq = z;
params.y_natural = y_natural;

fprintf('--- 初期条件 (自立モデル) ---\n');
fprintf('フレーム押し込み量: %.3f mm\n', x_push * 1000);
fprintf('計算された初期たわみ: %.3f mm\n', z * 1000);
fprintf('--------------------------------------\n');

% -------------------------------------------------------------------------
% 3. 周波数スイープの実行
% -------------------------------------------------------------------------
fprintf('周波数スイープを実行中...\n');

% 初期条件 (平衡点で静止した状態からスタート)
x0 = [z_target_initial; y_eq; 0; 0]; 
amp_data = zeros(size(freq_range));
h_wait = waitbar(0, '周波数スイープ中...');

for i = 1:length(freq_range)
    f_curr = freq_range(i);
    waitbar(i/length(freq_range), h_wait, sprintf('解析中: %.0f Hz', f_curr));
    
    params_sweep = params;
    params_sweep.freq = f_curr;
    
    % スイープ時は音波を最初からON (t_sound_on = 0)
    [~, x_sw] = ode45(@(t,x) equations_2DOF_SelfSustain(t, x, params_sweep, 0), tspan_sweep, x0);
    
    z_sw = x_sw(:,1);
    % 振動中心からの振幅を計算
    z_steady = z_sw(round(end/2):end); 
    amp_data(i) = (max(z_steady) - min(z_steady)) / 2;
end
close(h_wait);

[max_amp, idx_res] = max(amp_data);
f_res = freq_range(idx_res);

fprintf('スイープ完了。共振周波数: %.2f Hz\n', f_res);

% -------------------------------------------------------------------------
% 4. 本番シミュレーション (共振周波数にて)
% -------------------------------------------------------------------------
fprintf('共振周波数での時間応答を計算中...\n');

params.freq = f_res; 
options = odeset('RelTol',1e-6, 'AbsTol',1e-9);
[t, x] = ode45(@(t,x) equations_2DOF_SelfSustain(t, x, params, t_sound_on), tspan_main, x0, options);

z = x(:,1);
y = x(:,2);

% -------------------------------------------------------------------------
% 5. プロット
% -------------------------------------------------------------------------
% --- Figure 1: 周波数応答 ---
figure('Position', [50, 100, 600, 400], 'Color', 'w');
plot(freq_range, amp_data*1000, 'b.-', 'LineWidth', 1.5);
hold on;
plot(f_res, max_amp*1000, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
title('周波数応答 (剛性一定・自立モデル)');
xlabel('周波数 [Hz]'); ylabel('Z方向 振幅 [mm]');
grid on; legend('応答曲線', sprintf('共振点: %.1f Hz', f_res));

% --- Figure 2: 時間応答 ---
figure('Position', [700, 100, 800, 600], 'Color', 'w');

% (1) Z方向 (膜)
subplot(2,1,1);
plot(t, z*1000, 'b-');
hold on; 
xline(t_sound_on, 'r--', 'Label', 'Sound ON');
yline(z_target_initial*1000, 'k:', 'Label', '初期位置');
hold off;
ylabel('Z 変位 [mm]'); 
title(sprintf('膜のたわみ (駆動: %.2f Hz)', f_res)); 
grid on;

% (2) Y方向 (フレーム)
subplot(2,1,2);
% 拡張量としてプロット (y_eq - y)
y_extension = (params.y_eq - y) * 1000;
plot(t, y_extension, 'r-', 'LineWidth', 1.5);
hold on;
xline(t_sound_on, 'r--');
yline(0, 'g:', 'LineWidth', 1.5, 'Label', '初期位置');
hold off;
ylabel('拡張量 [mm]'); 
title('フレームの変位 (Y): 下がれば収縮、上がれば拡張'); 
grid on;


% =========================================================================
% 2自由度 運動方程式 (自立維持・剛性一定・音波あり)
% =========================================================================
function dx = equations_2DOF_SelfSustain(t, x, p, t_sound_on)
    z  = x(1); y  = x(2); dz = x(3); dy = x(4);
    
    % --- 1. 剛性は一定 (初期値のまま) ---
    k1_curr = p.k1;

    % --- 2. 復元力の計算 ---
    coupling_term = y - (p.Gamma * z^2);
    
    F_restore_z = k1_curr * z - p.K_couple * coupling_term * (2 * p.Gamma * z);
    F_restore_y = p.k2_base * (y - p.y_natural) + p.K_couple * coupling_term;
    
    % --- 3. 外力 (音波) ---
    if t >= t_sound_on
        % ランプを入れてショックを和らげる
        ramp = min(1.0, (t - t_sound_on)/0.1);
        F_sound = ramp * p.P0 * (p.L * p.w1) * cos(2 * pi * p.freq * t);
    else
        F_sound = 0;
    end
    
    % --- 4. 運動方程式 ---
    ddz = (F_sound - p.c1 * dz - F_restore_z) / p.m1;
    ddy = (0 - p.c2 * dy - F_restore_y) / p.m2;
    
    dx = [dz; dy; ddz; ddy];
end