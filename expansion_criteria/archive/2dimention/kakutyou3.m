% =========================================================================
% 音波で拡張するメタマテリアル（プリロード方式）
% 安定な2自由度モデル - 物理的に正しい平衡点の実装
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------

tspan = [0 8];
freq_range = 10:1:200;

% ★ 目標平衡点の設定（例）
z_target = 1.2e-3;      % [m] 目標Z方向変位（例: 1.2 mm）
y_target = -0.1e-3;     % [m] 目標Y方向変位（例: -0.1 mm、収縮状態）

params = struct( ...
    'L',       30e-3,       ...
    'L_side',  60e-3,       ...
    'L_up',    15e-3,       ...
    'P0',      2,           ...
    'E1',      197e9,       ...
    'w1',      10e-3,       ...
    't1',      0.02e-3,     ...
    'roh1',    7.93e3,      ...
    'delta1',  0.01,        ...
    'w2',      7.5e-3,      ...
    't2',      1e-3,        ...
    'S2',      126.87e-6,   ...
    'roh2',    1.24e3,      ...
    'k2',      177.035,     ...
    'delta2',  0.0278       ...
);

% -------------------------------------------------------------------------
% 2. パラメータの導出
% -------------------------------------------------------------------------

params = generate_params(params);

% -------------------------------------------------------------------------
% 3. プリロード力の計算（重要）
% -------------------------------------------------------------------------
% 目標平衡点 (z_target, y_target) で力が釣り合うように、
% Y方向のプリロード力 F_preload を計算する

% 幾何学的拘束: Y = -2Z²/L から導かれる関係
% 平衡点において、以下の力の釣り合いが成立する必要がある：
% Z方向: k1 * z_target + (非線形項) = 0
% Y方向: k2 * y_target + (膜からの連成項) = F_preload

% 分母の計算（平衡点での値）
denom_eq = params.L^2 - 2*z_target^2;

% Z方向の復元力（平衡点での値）
F_z_eq = params.k1 * z_target + (4*params.k2/denom_eq) * z_target^3;

% Y方向の復元力（平衡点での値、膜の張力による連成項を含む）
F_y_spring = params.k2 * y_target;  % フレームのばね力
F_y_from_membrane = (4*params.k1/params.L) * z_target^2;  % 膜からの力

% 平衡条件: F_preload = F_y_spring + F_y_from_membrane - F_z_coupling
% ここで、F_z_coupling はZ方向の力がY方向に伝達される項
F_z_coupling_to_y = (8*params.k2*z_target/denom_eq) * z_target^2;

% 必要なプリロード力（これを常に加えることで平衡点を維持）
F_preload = F_y_spring + F_y_from_membrane;

% パラメータ構造体に追加
params.z_eq = z_target;
params.y_eq = y_target;
params.F_preload = F_preload;

fprintf('========================================\n');
fprintf('目標平衡点:\n');
fprintf('  Z方向: %.4f mm\n', params.z_eq*1000);
fprintf('  Y方向: %.4f mm\n', params.y_eq*1000);
fprintf('  必要プリロード力: %.4f N\n', params.F_preload);
fprintf('========================================\n\n');

% -------------------------------------------------------------------------
% 4. 周波数スイープ
% -------------------------------------------------------------------------

fprintf('周波数スイープ実行中...\n');
amp = freq_sweep_2DOF_preload(params, freq_range, tspan);

[~, idx] = max(amp);
f_resonance = freq_range(idx);
omega_resonance = 2*pi*f_resonance;

fprintf('固有振動数: %.2f Hz\n\n', f_resonance);

% -------------------------------------------------------------------------
% 5. 時間応答シミュレーション
% -------------------------------------------------------------------------

fprintf('時間応答計算中...\n');

% 初期条件: 平衡点から静止状態でスタート
x0 = [params.z_eq; params.y_eq; 0; 0];

[t, x] = ode45(@(t,x) nonlinearForcedODE_2DOF_preload(t,x,params,omega_resonance), ...
               tspan, x0);

z = x(:,1);
y = x(:,2);

fprintf('計算完了\n\n');

% -------------------------------------------------------------------------
% 6. 結果のプロット
% -------------------------------------------------------------------------

% (1) 周波数応答
figure('Position', [100 100 800 600]);
plot(freq_range, amp*1000, 'b-', 'LineWidth', 2);
hold on;
plot(f_resonance, amp(idx)*1000, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('周波数 [Hz]');
ylabel('振幅 [mm]');
title(sprintf('周波数応答（固有振動数: %.2f Hz）', f_resonance));
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

% (2) Z方向の時間応答
figure('Position', [100 100 800 600]);
plot(t, z*1000, 'b-', 'LineWidth', 1.5);
hold on;
yline(params.z_eq*1000, 'r--', 'LineWidth', 2, 'Label', '平衡位置');
xlabel('時間 [s]');
ylabel('Z方向変位 [mm]');
title('Z方向（膜）の時間応答');
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
legend('z(t)', '平衡位置');

% (3) Y方向の時間応答
figure('Position', [100 100 800 600]);
plot(t, y*1000, 'r-', 'LineWidth', 1.5);
hold on;
yline(params.y_eq*1000, 'b--', 'LineWidth', 2, 'Label', '平衡位置');
xlabel('時間 [s]');
ylabel('Y方向変位 [mm]');
title('Y方向（フレーム）の時間応答');
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
legend('y(t)', '平衡位置');

% (4) 位相空間プロット
figure('Position', [100 100 800 600]);
plot(z*1000, y*1000, 'b-', 'LineWidth', 1);
hold on;
plot(params.z_eq*1000, params.y_eq*1000, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Z方向変位 [mm]');
ylabel('Y方向変位 [mm]');
title('位相空間プロット (Z-Y平面)');
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
legend('軌跡', '平衡点');

fprintf('全てのプロット完了\n');

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

function amp = freq_sweep_2DOF_preload(p, freq_range, tspan)
    amp = zeros(size(freq_range));
    x0 = [p.z_eq; p.y_eq; 0; 0];
    
    for i = 1:length(freq_range)
        omega = 2*pi*freq_range(i);
        [~, x] = ode45(@(t,x) nonlinearForcedODE_2DOF_preload(t,x,p,omega), ...
                       tspan, x0);
        
        z = x(:,1);
        z_deviation = z - p.z_eq;
        amp(i) = max(abs(z_deviation(end-100:end)));
    end
end

% -------------------------------------------------------------------------
% ★★★ プリロード方式の運動方程式（安定版）★★★
% -------------------------------------------------------------------------
function dx = nonlinearForcedODE_2DOF_preload(t, x, p, omega)
    % 状態変数
    z  = x(1);
    y  = x(2);
    dz = x(3);
    dy = x(4);
    
    % --- 外力 ---
    % Z方向の音波（正弦波加振）
    F_z_sound = p.P0 * (p.L - 2*p.t2) * p.w1 * cos(omega * t);
    
    % Y方向のプリロード力（常に一定、平衡点を維持するため）
    F_y_preload = p.F_preload;
    
    % --- 幾何学的非線形項の分母（特異点回避） ---
    denom = p.L^2 - 2*z^2;
    if abs(denom) < 1e-10
        denom = sign(denom) * 1e-10;
    end
    
    % --- 復元力（原点基準、元のモデルと同じ） ---
    
    % Z方向の復元力
    % (1) 膜の線形復元力
    F_z_linear = p.k1 * z;
    
    % (2) 幾何学的非線形連成項
    F_z_nonlinear = (4 * p.k2 / denom) * z^3;
    
    F_z_restore = F_z_linear + F_z_nonlinear;
    
    % Y方向の復元力
    % (1) フレームの線形復元力
    F_y_linear = p.k2 * y;
    
    % (2) 膜の張力による連成
    F_y_from_z = (4 * p.k1 / p.L) * z^2;
    
    F_y_restore = F_y_linear + F_y_from_z;
    
    % --- 見かけの質量（連成効果） ---
    m1_eff = p.m1 + (8 * p.m2 * z^2) / denom;
    m2_eff = p.m2 + (4 * p.m1 * z^2) / p.L;
    
    % --- 見かけの減衰 ---
    c1_eff = p.c1 + (8 * p.c2 * z^2) / denom;
    c2_eff = p.c2 + (4 * p.c1 * z^2) / p.L;
    
    % --- 慣性連成項 ---
    F_z_inertial = (8 * p.m2 * z / denom) * dz^2;
    F_y_inertial = (4 * p.m1 * z / p.L) * dz^2;
    
    % --- 加速度（プリロード力を含む） ---
    % Z方向の運動方程式（音波による加振）
    ddz = (F_z_sound - c1_eff*dz - F_z_restore - F_z_inertial) / m1_eff;
    
    % Y方向の運動方程式（プリロード力が復元力と釣り合う）
    ddy = (F_y_preload - c2_eff*dy - F_y_restore - F_y_inertial) / m2_eff;
    
    dx = [dz; dy; ddz; ddy];
end
