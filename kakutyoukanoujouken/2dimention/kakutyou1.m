% =========================================================================
% 音波で拡張するメタマテリアルの2自由度モデル。z方向が０中心になっている。
% 2-DOF Model for Acoustic Wave-Expanded Metamaterial
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------

tspan = [0 8];              % シミュレーション時間 [s]
freq_range = 10:1:200;      % 周波数スイープ範囲 [Hz]

% 初期たわみ設定（拡張モデルの特徴）
y_offset = -0.1e-3;         % 初期たわみ位置 [m] (負の値 = 収縮状態)

% 物理パラメータ
params = struct( ...
    'L',       30e-3,       ... % 膜の長さ、枠の縦長さ [m]
    'L_side',  60e-3,       ... % 枠の横長さ [m]
    'L_up',    15e-3,       ... % 上底の長さ [m]
    'P0',      2,           ... % 外力振幅（音圧） [Pa]
    'E1',      197e9,       ... % 膜のヤング率 [Pa] (SUS304)
    'w1',      10e-3,       ... % 膜の幅 [m]
    't1',      0.02e-3,     ... % 膜の厚さ [m]
    'roh1',    7.93e3,      ... % 膜の密度 [kg/m^3] (SUS304)
    'delta1',  0.01,        ... % 膜の対数減衰率
    'w2',      7.5e-3,      ... % 枠の幅 [m]
    't2',      1e-3,        ... % 枠の厚さ [m]
    'S2',      126.87e-6,   ... % 枠の断面積 [m^2]
    'roh2',    1.24e3,      ... % 枠の密度 [kg/m^3]
    'k2',      177.035,     ... % 枠のばね定数 [N/m] (ANSYS解析値)
    'delta2',  0.0278,      ... % 枠の対数減衰率
    'y_offset', y_offset    ... % 初期たわみ位置
);

% -------------------------------------------------------------------------
% 2. パラメータの導出
% -------------------------------------------------------------------------

params = generate_params(params);

% 初期たわみ状態での平衡変位を計算
z_equilibrium = sqrt(-params.y_offset * params.L / 2);

fprintf('========================================\n');
fprintf('初期たわみ設定\n');
fprintf('  Y方向オフセット: %.4f mm\n', params.y_offset*1000);
fprintf('  対応するZ方向変位: %.4f mm\n', z_equilibrium*1000);
fprintf('========================================\n\n');

% -------------------------------------------------------------------------
% 3. 周波数スイープ（固有振動数の探索）
% -------------------------------------------------------------------------

fprintf('周波数スイープ実行中...\n');
amp = freq_sweep_2DOF(params, freq_range, tspan, z_equilibrium);

[~, idx] = max(amp);
f_resonance = freq_range(idx);
omega_resonance = 2*pi*f_resonance;

fprintf('固有振動数: %.2f Hz\n\n', f_resonance);

% -------------------------------------------------------------------------
% 4. 共振周波数での時間応答シミュレーション
% -------------------------------------------------------------------------

fprintf('共振周波数での時間応答計算中...\n');

% 初期条件: 初期たわみ状態（静止）
x0 = [z_equilibrium; params.y_offset; 0; 0]; % [z; y; dz; dy]

% ODE求解
[t, x] = ode45(@(t,x) nonlinearForcedODE_2DOF(t,x,params,omega_resonance), ...
               tspan, x0);

z = x(:,1);  % Z方向変位
y = x(:,2);  % Y方向変位
dz = x(:,3); % Z方向速度
dy = x(:,4); % Y方向速度

fprintf('計算完了\n\n');

% -------------------------------------------------------------------------
% 5. 結果のプロット
% -------------------------------------------------------------------------

% (1) 周波数応答
figure('Position', [100 100 800 600]);
plot(freq_range, amp*1000, 'b-', 'LineWidth', 2);
hold on;
plot(f_resonance, amp(idx)*1000, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('周波数 [Hz]');
ylabel('振幅 [mm]');
title(sprintf('周波数応答（固有振動数: %.2f Hz）', f_resonance));
legend('振幅', sprintf('共振点: %.2f Hz', f_resonance));
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% (2) Z方向の時間応答
figure('Position', [100 100 800 600]);
plot(t, z*1000, 'b-', 'LineWidth', 1.5);
hold on;
yline(z_equilibrium*1000, 'r--', 'LineWidth', 2, 'Label', '初期平衡位置');
xlabel('時間 [s]');
ylabel('Z方向変位 [mm]');
title('Z方向（膜）の時間応答');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('z(t)', '初期平衡位置');

% (3) Y方向の時間応答
figure('Position', [100 100 800 600]);
plot(t, y*1000, 'r-', 'LineWidth', 1.5);
hold on;
yline(params.y_offset*1000, 'b--', 'LineWidth', 2, 'Label', '初期たわみ位置');
xlabel('時間 [s]');
ylabel('Y方向変位 [mm]');
title('Y方向（フレーム）の時間応答');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('y(t)', '初期たわみ位置');

% (4) 位相空間プロット (Z-Y平面)
figure('Position', [100 100 800 600]);
plot(z*1000, y*1000, 'b-', 'LineWidth', 1);
hold on;
plot(z_equilibrium*1000, params.y_offset*1000, 'ro', ...
     'MarkerSize', 10, 'LineWidth', 2);
xlabel('Z方向変位 [mm]');
ylabel('Y方向変位 [mm]');
title('位相空間プロット (Z-Y平面)');
legend('軌跡', '初期平衡点');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
grid on;
axis equal;

fprintf('全てのプロット完了\n');

% =========================================================================
% 関数定義
% =========================================================================

% -------------------------------------------------------------------------
% パラメータ導出関数
% -------------------------------------------------------------------------
function p = generate_params(p)
    % 膜のパラメータ
    p.I1 = p.w1 * p.t1^3 / 12;              % 断面二次モーメント [m^4]
    p.k1 = 384 * p.E1 * p.I1 / (p.L^3);     % ばね定数 [N/m]
    p.m1 = p.roh1 * p.L * p.w1 * p.t1;      % 質量 [kg]
    p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1); % 減衰係数 [Ns/m]
    
    % 枠のパラメータ
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3); % 断面積 [m^2]
    p.m2 = p.S2 * p.w2 * p.roh2;            % 質量 [kg]
    p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2); % 減衰係数 [Ns/m]
    
    fprintf('導出パラメータ:\n');
    fprintf('  膜: m1=%.4e kg, k1=%.2f N/m, c1=%.4e Ns/m\n', p.m1, p.k1, p.c1);
    fprintf('  枠: m2=%.4e kg, k2=%.2f N/m, c2=%.4e Ns/m\n', p.m2, p.k2, p.c2);
end

% -------------------------------------------------------------------------
% 周波数スイープ関数（2自由度版）
% -------------------------------------------------------------------------
function amp = freq_sweep_2DOF(p, freq_range, tspan, z_eq)
    amp = zeros(size(freq_range));
    
    % 初期たわみ状態から開始
    x0 = [z_eq; p.y_offset; 0; 0];
    
    for i = 1:length(freq_range)
        omega = 2*pi*freq_range(i);
        
        % ODE求解
        [~, x] = ode45(@(t,x) nonlinearForcedODE_2DOF(t,x,p,omega), ...
                       tspan, x0);
        
        % Z方向変位の振幅（平衡点からの偏差）
        z = x(:,1);
        z_deviation = z - z_eq;
        amp(i) = max(abs(z_deviation(end-100:end)));
    end
end

% -------------------------------------------------------------------------
% 2自由度連立運動方程式
% -------------------------------------------------------------------------
function dx = nonlinearForcedODE_2DOF(t, x, p, omega)
    % 状態変数の展開
    z  = x(1);  % Z方向変位（膜）
    y  = x(2);  % Y方向変位（フレーム）
    dz = x(3);  % Z方向速度
    dy = x(4);  % Y方向速度
    
    % --- 外力定義 ---
    % Z方向の音波（正弦波加振）
    F_z_sound = p.P0 * (p.L - 2*p.t2) * p.w1 * cos(omega * t);
    
    % Y方向の外力（通常は0、初期たわみは復元力で表現）
    F_y_external = 0;
    
    % --- 幾何学的非線形項の分母（特異点回避） ---
    denom = p.L^2 - 2*z^2;
    if abs(denom) < 1e-10
        denom = sign(denom) * 1e-10;
    end
    
    % --- Z方向の復元力 ---
    % (1) 膜自身の線形復元力
    F_z_linear = p.k1 * z;
    
    % (2) フレームとの連成による非線形復元力
    %     幾何学的拘束 Y = -2Z²/L から導出
    F_z_coupling = (4 * p.k2 / denom) * z^3;
    
    % 合計復元力
    F_z_restore = F_z_linear + F_z_coupling;
    
    % --- Y方向の復元力 ---
    % (1) フレームの線形復元力（初期たわみ位置からの変位）
    F_y_linear = p.k2 * (y - p.y_offset);
    
    % (2) 膜の張力による連成力
    %     膜の変位がY方向に与える影響
    F_y_coupling = (4 * p.k1 / p.L) * z^2 * sign(y);
    
    % 合計復元力
    F_y_restore = F_y_linear + F_y_coupling;
    
    % --- 見かけの質量（連成効果） ---
    m1_eff = p.m1 + (8 * p.m2 * z^2) / denom;
    m2_eff = p.m2 + (4 * p.m1 * z^2) / p.L;
    
    % --- 見かけの減衰（連成効果） ---
    c1_eff = p.c1 + (8 * p.c2 * z^2) / denom;
    c2_eff = p.c2 + (4 * p.c1 * z^2) / p.L;
    
    % --- 慣性連成項（速度の2乗項） ---
    % Coriolis/centrifugal-like terms from geometric coupling
    F_z_inertial = (8 * p.m2 * z / denom) * dz^2;
    F_y_inertial = (4 * p.m1 * z / p.L) * dz^2;
    
    % --- 加速度の計算 ---
    % Z方向の運動方程式
    ddz = (F_z_sound - c1_eff*dz - F_z_restore - F_z_inertial) / m1_eff;
    
    % Y方向の運動方程式
    ddy = (F_y_external - c2_eff*dy - F_y_restore - F_y_inertial) / m2_eff;
    
    % 状態ベクトルの微分
    dx = [dz; dy; ddz; ddy];
end
