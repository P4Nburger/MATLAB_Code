% =========================================================================
% untitled.m の2自由度版
% Z方向とY方向を独立した自由度として扱う
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------

% --- 各荷重のタイミングと ramp 時間を設定 ---
t_ramp_static = 1;          % [s] 静的荷重を最大値まで増加させる時間
t_start_dynamic = 2;        % [s] 音波（動的荷重）をかけ始める時間
t_ramp_dynamic = 1.0;       % [s] 動的荷重を1秒かけて最大振幅にする

% --- 目標とする静的変位（初期たわみ）をここで設定 ---
y_static_target = -0.05e-3; % [m] Y方向の目標変位（例: -0.05 mm）

tspan = [0 8];
freq_range = 10:1:200;

params = struct( ...
    'L',       30e-3,       ...
    'L_side',  60e-3,       ...
    'P0',      2,           ...
    'E1',      197e9,       ...
    'w1',      10e-3,       ...
    't1',      0.02e-3,     ...
    'roh1',    7.93e3,      ...
    'delta1',  0.01,        ...
    'w2',      7.5e-3,      ...
    'L_up',    15e-3,       ...
    't2',      1e-3,        ...
    'S2',      126.87e-6,   ...
    'roh2',    1.24e3,      ...
    'k2',      177.035,     ...
    'delta2',  0.0278       ...
);

% -------------------------------------------------------------------------
% 2. 荷重計算とパラメータ設定
% -------------------------------------------------------------------------

params = generate_params(params);

if y_static_target > 0
    error('y_static_targetは物理的に負の値である必要があります。');
end

% Z方向の目標変位を計算（幾何学的拘束 Y = -2Z²/L から）
z_static_target = sqrt(-y_static_target * params.L / 2);

% 目標平衡点での必要な静的荷重を計算
z_eq = z_static_target;
y_eq = y_static_target;
denom_eq = params.L^2 - 2*z_eq^2;

% Z方向の静的荷重（膜を押すため）
F_z_static_required = params.k1*z_eq + (4*params.k2/denom_eq)*z_eq^3;

% Y方向の静的荷重（フレームを引っ張るため）
F_y_static_required = params.k2*y_eq + (4*params.k1/params.L)*z_eq^2;

% パラメータに追加
params.F_z_static_max = F_z_static_required;
params.F_y_static_max = F_y_static_required;
params.t_ramp_static = t_ramp_static;
params.t_start_dynamic = t_start_dynamic;
params.t_ramp_dynamic = t_ramp_dynamic;
params.z_static_target = z_static_target;
params.y_static_target = y_static_target;

% 初期条件: [z; y; dz; dy] = [0; 0; 0; 0] からスタート
x0 = [0; 0; 0; 0];

fprintf('---------------------------------------------------\n');
fprintf('目標変位:\n');
fprintf('  Z方向: %.4f mm\n', z_static_target*1000);
fprintf('  Y方向: %.4f mm\n', y_static_target*1000);
fprintf('必要な静的荷重:\n');
fprintf('  Z方向: %.4f N\n', F_z_static_required);
fprintf('  Y方向: %.4f N\n', F_y_static_required);
fprintf('タイミング:\n');
fprintf('  静的荷重 ramp時間: %.1f 秒\n', t_ramp_static);
fprintf('  動的荷重 開始時間: %.1f 秒\n', t_start_dynamic);
fprintf('  動的荷重 ramp時間: %.1f 秒\n', t_ramp_dynamic);
fprintf('---------------------------------------------------\n');

% -------------------------------------------------------------------------
% 3. 周波数スイープ
% -------------------------------------------------------------------------

fprintf('周波数スイープ実行中...\n');
amp = freq_sweep_2DOF(params, freq_range, tspan, x0);

[~, idx] = max(amp);
f_res = freq_range(idx);
omega_res = 2*pi*f_res;

fprintf('シミュレーション結果:\n');
fprintf('  固有振動数（目標変位周り）: %.2f Hz\n', f_res);

% -------------------------------------------------------------------------
% 4. 時間応答計算
% -------------------------------------------------------------------------

fprintf('時間応答計算中...\n');
[t, x] = ode45(@(t,x) nonlinearForcedODE_2DOF(t,x,params,omega_res), tspan, x0);

z = x(:,1);
y = x(:,2);
dz = x(:,3);
dy = x(:,4);

fprintf('計算完了\n');

% -------------------------------------------------------------------------
% 5. プロット
% -------------------------------------------------------------------------

% (1) 周波数応答
figure('Position', [100 100 800 600]);
plot(freq_range, amp*1000, 'b', 'LineWidth', 2);
title('振幅応答（静的変位周り）');
xlabel('周波数 [Hz]');
ylabel('振幅 [mm]');
set(gca, 'Fontsize', 20, 'FontWeight', 'bold');
grid on;

% (2) Z方向の時間応答
figure('Position', [100 100 800 600]);
plot(t, z*1000, 'b-', 'LineWidth', 1.5);
hold on;
yline(z_static_target * 1000, 'r--', 'LineWidth', 2, 'Label', '目標平衡位置');
hold off;
title('時間応答：Z方向（膜）');
xlabel('時間 [s]');
ylabel('変位 z(t) [mm]');
set(gca, 'Fontsize', 20, 'FontWeight', 'bold');
grid on;

% (3) Y方向の時間応答
figure('Position', [100 100 800 600]);
plot(t, y*1000, 'r-', 'LineWidth', 1.5);
hold on;
yline(y_static_target * 1000, 'b--', 'LineWidth', 2, 'Label', '目標平衡位置');
hold off;
title('時間応答：Y方向（フレーム）');
xlabel('時間 [s]');
ylabel('変位 y(t) [mm]');
set(gca, 'Fontsize', 20, 'FontWeight', 'bold');
grid on;

% (4) 位相空間プロット (Z-Y平面)
figure('Position', [100 100 800 600]);
plot(z*1000, y*1000, 'b-', 'LineWidth', 1);
hold on;
plot(z_static_target*1000, y_static_target*1000, 'ro', ...
     'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '目標平衡点');
xlabel('Z方向変位 [mm]');
ylabel('Y方向変位 [mm]');
title('位相空間プロット (Z-Y平面)');
set(gca, 'Fontsize', 20, 'FontWeight', 'bold');
grid on;
legend;

% (5) 印加荷重の時間変化
F_z_static_history = params.F_z_static_max * min(1.0, t / params.t_ramp_static);
F_y_static_history = params.F_y_static_max * min(1.0, t / params.t_ramp_static);

t_since_dynamic_start = t - params.t_start_dynamic;
ramp_dynamic_factor = min(1.0, max(0, t_since_dynamic_start / params.t_ramp_dynamic));
F_dynamic_envelope = ramp_dynamic_factor * params.P0 * (params.L - 2*params.t2) * params.w1;

figure('Position', [100 100 800 600]);
hold on;
plot(t, F_z_static_history, 'm-', 'LineWidth', 2.5, 'DisplayName', 'Z方向静的荷重');
plot(t, F_y_static_history, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Y方向静的荷重');
plot(t, F_dynamic_envelope, 'c-', 'LineWidth', 2.5, 'DisplayName', '動的荷重の振幅');
hold off;
title('印加される荷重の時間変化');
xlabel('時間 [s]');
ylabel('荷重 [N]');
legend('Location', 'best');
set(gca, 'Fontsize', 20, 'FontWeight', 'bold');
grid on;

% =========================================================================
% 関数定義
% =========================================================================

function p = generate_params(p)
    p.I1 = p.w1 * p.t1^3 / 12;
    p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.m1 = p.roh1 * p.L * p.w1 * p.t1;
    p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1);
    
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3);
    p.m2 = p.S2 * p.w2 * p.roh2;
    p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2);
    
    fprintf('導出パラメータ:\n');
    fprintf('  膜: m1=%.4e kg, k1=%.2f N/m, c1=%.4e Ns/m\n', p.m1, p.k1, p.c1);
    fprintf('  枠: m2=%.4e kg, k2=%.2f N/m, c2=%.4e Ns/m\n', p.m2, p.k2, p.c2);
end

function amp = freq_sweep_2DOF(p, freq_range, tspan, x0)
    amp = zeros(size(freq_range));
    
    for i = 1:length(freq_range)
        omega = 2*pi*freq_range(i);
        [~, x] = ode45(@(t,x) nonlinearForcedODE_2DOF(t,x,p,omega), tspan, x0);
        
        z = x(:,1);
        z_offset = z - p.z_static_target;
        amp(i) = max(abs(z_offset(end-100:end)));
    end
end

% -------------------------------------------------------------------------
% ★ 2自由度運動方程式 ★
% -------------------------------------------------------------------------
function dx = nonlinearForcedODE_2DOF(t, x, p, omega)
    % 状態変数の展開
    z  = x(1);  % Z方向変位
    y  = x(2);  % Y方向変位
    dz = x(3);  % Z方向速度
    dy = x(4);  % Y方向速度
    
    % --- 静的荷重（ramp関数で徐々に増加） ---
    ramp_static_factor = min(1.0, t / p.t_ramp_static);
    F_z_static_current = p.F_z_static_max * ramp_static_factor;
    F_y_static_current = p.F_y_static_max * ramp_static_factor;
    
    % --- 動的荷重（Z方向のみ、指定時刻以降に加振） ---
    F_z_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2*p.t2) * p.w1;
        F_z_dynamic = amplitude * cos(omega * time_since_start);
    end
    
    % Z方向の総外力
    F_z_total = F_z_static_current + F_z_dynamic;
    
    % Y方向の総外力
    F_y_total = F_y_static_current;
    
    % --- 幾何学的非線形項の分母（特異点回避） ---
    denom = p.L^2 - 2*z^2;
    if abs(denom) < 1e-10
        denom = sign(denom) * 1e-10;
    end
    
    % --- 復元力 ---
    % Z方向の復元力
    F_z_restore = p.k1 * z + (4*p.k2/denom) * z^3;
    
    % Y方向の復元力
    F_y_restore = p.k2 * y + (4*p.k1/p.L) * z^2;
    
    % --- 見かけの質量（連成効果） ---
    m1_eff = p.m1 + (8*p.m2*z^2) / denom;
    m2_eff = p.m2 + (4*p.m1*z^2) / p.L;
    
    % --- 見かけの減衰 ---
    c1_eff = p.c1 + (8*p.c2*z^2) / denom;
    c2_eff = p.c2 + (4*p.c1*z^2) / p.L;
    
    % --- 慣性連成項 ---
    F_z_inertial = (8*p.m2*z/denom) * dz^2;
    F_y_inertial = (4*p.m1*z/p.L) * dz^2;
    
    % --- 加速度の計算 ---
    ddz = (F_z_total - c1_eff*dz - F_z_restore - F_z_inertial) / m1_eff;
    ddy = (F_y_total - c2_eff*dy - F_y_restore - F_y_inertial) / m2_eff;
    
    % 状態ベクトルの微分
    dx = [dz; dy; ddz; ddy];
end
