% =========================================================================
% 2自由度 シミュレーション (自立維持 + 剛性一定 + 周波数スイープ)
% ・初期たわみ: フレームの内力(y_natural)で自立維持（真の平衡点）
% ・剛性変化: ★なし (k1は常に一定)★
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
% --- スイープ設定 ---
freq_range   = 1:1:150;   % [Hz]
tspan_sweep  = [0 2.0];   % [s]

% --- 音波設定 ---
P0_drive     = 2;         % [Pa]

% --- 本番シミュレーション設定 ---
tspan_main   = [0 5.0];   % [s]
t_sound_on   = 1;       % [s]

% --- 初期たわみ設定 ---0.02でスナップスルー,0.03でただの振動(12/4時点)
x_push = 0.03e-3;          % [m] フレーム押し込み量
C_eff  = 0.2794;          % 有効形状定数

% x_push = (1/2)*C_eff*delta_0^2  [mm系] を m系に直した実装
z_target_initial_mm = sqrt(2 * (x_push * 1000) / C_eff);  % [mm]
z_target_initial    = z_target_initial_mm / 1000;         % [m]

% --- 構造・材料パラメータ ---
params = struct( ...
    'L',   18e-3, ...      % [m]
    'w1',  5e-3, ...      % [m]
    't1',  0.1e-3, ...    % [m]
    'E1',  3.45e9, ...      % [Pa]
    'rho1',1250, ...       % [kg/m^3]
    'S2',  152.522e-6, ...  % [m^2]
    'rho2',1250, ...       % [kg/m^3]
    'k2_base',222.0, ...   % [N/m]
    'delta',0.05, ...      % 減衰比
    'P0',  P0_drive, ...   % [Pa]
    'freq',0 ...           % [Hz]（後で上書き）
);

% -------------------------------------------------------------------------
% 2. 物理定数と「真の平衡点」の計算
% -------------------------------------------------------------------------
params.m1 = params.rho1 * params.L * params.w1 * params.t1;
params.m2 = params.rho2 * params.S2 * params.L;
params.I1 = params.w1 * params.t1^3 / 12;

params.k1       = 384 * params.E1 * params.I1 / (params.L^3);
params.K_couple = params.E1 * (params.w1 * params.t1) / params.L;
params.c1       = 2 * params.delta * sqrt(params.m1 * params.k1);
params.c2       = 2 * params.delta * sqrt(params.m2 * params.k2_base);

% --- 幾何学係数 Gamma ---
% x_push = (1/2)*C_eff*delta_0^2  と  y ≈ Gamma*z^2 の対応より Gamma ≈ (C_eff/2)×(mm→m換算)
params.Gamma = (C_eff / 2) * 1000;   % [1/m]

% --- 欲しい初期たわみを z_eq として固定 ---
z_eq = z_target_initial;

% ラグランジュからの厳密な釣り合い条件：
% Z方向: 0 = k1*z_eq - 2*Gamma*z_eq*K_couple*(y_eq - Gamma*z_eq^2)
%  → coupling_eq = y_eq - Gamma*z_eq^2 = k1 / (2*Gamma*K_couple)
coupling_eq = params.k1 / (2 * params.Gamma * params.K_couple);

% したがって平衡点の y_eq は
y_eq = params.Gamma * z_eq^2 + coupling_eq;

% Y方向: 0 = k2*(y_eq - y_natural) + K_couple*(y_eq - Gamma*z_eq^2)
%  → y_natural = y_eq + (K_couple/k2)*coupling_eq
y_natural = y_eq + (params.K_couple * coupling_eq) / params.k2_base;

params.z_eq      = z_eq;
params.y_eq      = y_eq;
params.y_natural = y_natural;

fprintf('--- 初期条件 (自立モデル: 真の平衡点) ---\n');
fprintf('フレーム押し込み量:    %.3f mm\n', x_push * 1000);
fprintf('目標初期たわみ z_eq:    %.3f mm\n', z_eq * 1000);
fprintf('対応する平衡点 y_eq:    %.3f mm\n', y_eq * 1000);
fprintf('フレーム自然長 y_nat:   %.3f mm\n', y_natural * 1000);
fprintf('--------------------------------------\n');

% 初期条件（真の力の釣り合い点で静止からスタート）
x0 = [z_eq; y_eq; 0; 0];

% -------------------------------------------------------------------------
% 3. 周波数スイープの実行
% -------------------------------------------------------------------------
fprintf('周波数スイープを実行中...\n');

amp_data = zeros(size(freq_range));
h_wait   = waitbar(0, '周波数スイープ中...');

for i = 1:length(freq_range)
    f_curr = freq_range(i);
    waitbar(i/length(freq_range), h_wait, sprintf('解析中: %.0f Hz', f_curr));
    
    params_sweep      = params;
    params_sweep.freq = f_curr;
    
    % スイープ時は音波を最初からON (t_sound_on = 0)
    [~, x_sw] = ode45(@(t,x) equations_2DOF_SelfSustain(t, x, params_sweep, 0), ...
                      tspan_sweep, x0);
    
    z_sw    = x_sw(:,1);
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
options     = odeset('RelTol',1e-6, 'AbsTol',1e-9);

[t, x] = ode45(@(t,x) equations_2DOF_SelfSustain(t, x, params, t_sound_on), ...
               tspan_main, x0, options);

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
yline(z_eq*1000, 'k:', 'Label', '初期位置');
hold off;
ylabel('z(t)[mm]'); 
title(sprintf('膜のたわみ (Z)(駆動: %.2f Hz)', f_res)); 
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
ylabel('y(t)[mm]'); 
title('フレームの変位 (Y)'); 
grid on;

% =========================================================================
% 2自由度 運動方程式 (自立維持・剛性一定・音波あり)
% =========================================================================
function dx = equations_2DOF_SelfSustain(t, x, p, t_sound_on)
    z  = x(1); 
    y  = x(2); 
    dz = x(3); 
    dy = x(4);
    
    % --- 1. 剛性は一定 (初期値のまま) ---
    k1_curr = p.k1;

    % --- 2. 復元力の計算 ---
    coupling_term = y - (p.Gamma * z^2);
    
    F_restore_z = k1_curr * z ...
                  - p.K_couple * coupling_term * (2 * p.Gamma * z);
    F_restore_y = p.k2_base * (y - p.y_natural) ...
                  + p.K_couple * coupling_term;
    
    % --- 3. 外力 (音波) ---
    if t >= t_sound_on
        % ランプを入れてショックを和らげる
        ramp    = min(1.0, (t - t_sound_on)/0.1);
        F_sound = ramp * p.P0 * (p.L * p.w1) * cos(2 * pi * p.freq * t);
    else
        F_sound = 0;
    end
    
    % --- 4. 運動方程式 ---
    ddz = (F_sound - p.c1 * dz - F_restore_z) / p.m1;
    ddy = (0 - p.c2 * dy - F_restore_y)     / p.m2;
    
    dx = [dz; dy; ddz; ddy];
end
