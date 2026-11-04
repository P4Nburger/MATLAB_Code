% =========================================================================
% 2自由度 シミュレーション (v2改ベース: 剛性変化なし・純粋共振モデル)
% ・初期たわみ: フレーム押し込み量から算出
% ・周波数スイープ: あり
% ・剛性変化: ★無効化 (k1は常に一定)★
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
% --- スイープ設定 ---
freq_range = 10:1:200;  % [Hz] スイープする周波数範囲
tspan_sweep = [0 1.0];  % [s] スイープ時の計算時間

% --- 本番シミュレーション設定 ---
tspan_main = [0 1.0];   % [s] 時間応答の計算時間
t_sound_on = 0.2;       % [s] 音波開始

% --- 初期たわみ設定 (実験値準拠) ---
x_push = 0.1e-3;        % [m] フレーム押し込み量
C_eff = 0.2794;         % 有効形状定数
z_target_initial_mm = sqrt(2 * (x_push * 1000) / C_eff);
z_target_initial = z_target_initial_mm / 1000; % [m]

% --- 音波設定 ---
P0_drive = 25;          % [Pa]

% --- 剛性変化パラメータ (今回は無効化されますが変数は残します) ---
stiffness_factor = 0;   % ★ 0に設定して無効化
stiffness_threshold = 0;
amplitude_sensitivity = 0;
amplitude_decay = 0;

params = struct( ...
    'L', 30e-3, 'w1', 10e-3, 't1', 0.02e-3, 'E1', 197e9, 'rho1', 7930, ...
    'S2', 126.87e-6, 'rho2', 1240, 'k2_base', 177.0, ...
    'delta', 0.05, 'P0', P0_drive, 'freq', 0, ... 
    'stiffness_factor', stiffness_factor, ...
    'threshold', stiffness_threshold, ...
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

% 幾何学係数 Gamma
params.Gamma = (C_eff / 2) * 1000; 

% --- 平衡点の計算 ---
z = z_target_initial;
y_eq = x_push; 

% 平衡点を作るための自然長 y_natural の逆算
geom_error = y_eq - (params.Gamma * z^2); 
F_tension_y_equivalent = params.k2_base * 0.1e-3; 
y_natural = y_eq + F_tension_y_equivalent / params.k2_base;

params.y_eq = y_eq;
params.z_eq = z;
params.y_natural = y_natural;

fprintf('--- 初期条件 ---\n');
fprintf('フレーム押し込み量: %.3f mm\n', x_push * 1000);
fprintf('計算された初期たわみ: %.3f mm\n', z * 1000);
fprintf('--------------------------------------\n');

% -------------------------------------------------------------------------
% 3. 周波数スイープの実行
% -------------------------------------------------------------------------
fprintf('周波数スイープを実行中 (%.0f Hz ～ %.0f Hz)...\n', min(freq_range), max(freq_range));

x0 = [z_target_initial; y_eq; 0; 0; 0]; % 状態変数: [z, y, dz, dy, A]
amp_data = zeros(size(freq_range));
h_wait = waitbar(0, '周波数スイープ中...');

for i = 1:length(freq_range)
    f_curr = freq_range(i);
    waitbar(i/length(freq_range), h_wait, sprintf('解析中: %.0f Hz', f_curr));
    
    params_sweep = params;
    params_sweep.freq = f_curr;
    
    % スイープ時は音波を最初からON (t_sound_on = 0)
    [~, x_sw] = ode45(@(t,x) equations_2DOF_Constant(t, x, params_sweep, 0), tspan_sweep, x0);
    
    z_sw = x_sw(:,1);
    z_steady = z_sw(round(end/2):end); 
    amp_data(i) = (max(z_steady) - min(z_steady)) / 2;
end
close(h_wait);

[max_amp, idx_res] = max(amp_data);
f_res = freq_range(idx_res);

fprintf('スイープ完了。\n');
fprintf('検出された共振周波数: %.2f Hz (最大振幅: %.3f mm)\n', f_res, max_amp*1000);
fprintf('--------------------------------------\n');

% -------------------------------------------------------------------------
% 4. 本番シミュレーション (共振周波数にて)
% -------------------------------------------------------------------------
fprintf('共振周波数での時間応答を計算中...\n');

params.freq = f_res; 
options = odeset('RelTol',1e-6, 'AbsTol',1e-9);
[t, x] = ode45(@(t,x) equations_2DOF_Constant(t, x, params, t_sound_on), tspan_main, x0, options);

z = x(:,1);
y = x(:,2);
A = x(:,5);

% -------------------------------------------------------------------------
% 5. プロット
% -------------------------------------------------------------------------
% --- Figure 1: 周波数応答 ---
figure('Position', [50, 100, 600, 400], 'Color', 'w');
plot(freq_range, amp_data*1000, 'b.-', 'LineWidth', 1.5);
hold on;
plot(f_res, max_amp*1000, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
title('周波数応答 (剛性変化なし)');
xlabel('周波数 [Hz]'); ylabel('Z方向 振幅 [mm]');
grid on; legend('応答曲線', sprintf('共振点: %.1f Hz', f_res));

% --- Figure 2: 時間応答 ---
figure('Position', [700, 100, 800, 900], 'Color', 'w');

% (1) Z方向 (膜)
subplot(3,1,1);
plot(t, z*1000, 'b-');
hold on; 
xline(t_sound_on, 'r--', 'Label', 'Sound ON');
yline(z_target_initial*1000, 'k:', 'Label', '初期位置');
hold off;
ylabel('Z 変位 [mm]'); 
title(sprintf('膜のたわみ (駆動: %.2f Hz)', f_res)); 
grid on;

% (2) Y方向 (フレーム)
subplot(3,1,2);
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

% (3) 内部変数 (剛性は一定のはず)
subplot(3,1,3);
yyaxis left; 
plot(t, A, 'b-', 'LineWidth', 1.5);
ylabel('振動の激しさ A');

yyaxis right; 
% 剛性は一定
k1_constant = ones(size(t)) * params.k1; 
plot(t, k1_constant, 'r-', 'LineWidth', 1.5);
ylabel('膜の剛性 k1 [N/m]');
ylim([params.k1*0.9, params.k1*1.1]);
title('剛性は一定 (変化なし)'); 
grid on; xlabel('時間 [s]');


% =========================================================================
% 2自由度 運動方程式 (剛性変化なし・一定剛性)
% =========================================================================
function dx = equations_2DOF_Constant(t, x, p, t_sound_on)
    z  = x(1);
    y  = x(2);
    dz = x(3);
    dy = x(4);
    A  = x(5);
    
    % --- 1. 剛性は一定 (変化なし) ---
    % ロジックを無効化し、初期値 k1 をそのまま使う
    k1_curr = p.k1; 

    % --- 2. 復元力の計算 (連成) ---
    coupling_term = y - (p.Gamma * z^2);
    
    % Z方向
    F_restore_z = k1_curr * z - p.K_couple * coupling_term * (2 * p.Gamma * z);
    
    % Y方向
    F_restore_y = p.k2_base * (y - p.y_natural) + p.K_couple * coupling_term;
    
    % --- 3. 外力 ---
    if t >= t_sound_on
        ramp = min(1.0, (t - t_sound_on)/0.1);
        F_sound = ramp * p.P0 * (p.L * p.w1) * cos(2 * pi * p.freq * t);
    else
        F_sound = 0;
    end
    
    % --- 4. 運動方程式 ---
    ddz = (F_sound - p.c1 * dz - F_restore_z) / p.m1;
    ddy = (0 - p.c2 * dy - F_restore_y) / p.m2;
    
    % --- 5. 振幅Aの更新 (計算だけはしておく) ---
    dA = p.amp_sens * abs(dz) - p.amp_decay * A;
    
    dx = [dz; dy; ddz; ddy; dA];
end