% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ 比較実験の条件をここで設定 ★ ---
% ケースBでかけるプリテンションの強さを定義
z_preload_for_tension = 0.2e-3; % [m] プリテンションをかけるための仮想的なZ変位

% 両ケース共通でかける音波の条件
f_drive = 120;            % [Hz] ★ 比較のための駆動周波数を固定 (共振点から少しずらすと差が分かりやすい)
P0_drive = 10;            % [Pa] ★ 音圧
t_start_dynamic = 1.0;    % [s] 音波をかけ"始める"時間
t_ramp_dynamic = 1.0;     % [s] 動的荷重をこの時間かけて最大振幅にする
% ----------------------------------------------------

tspan = [0 8];
params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', P0_drive, 'E1', 197e9, ...
        'w1', 10e-3, 't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.01, ...
        'w2', 7.5e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, ...
        'roh2', 1.24e3, 'k2', 177.035, 'delta2', 0.0278 ...
    );

% -------------------------
% 2. シミュレーション実行
% -------------------------
params_base = generate_params(params);
omega_drive = 2 * pi * f_drive;
x0 = [0;0];

% --- ケースA: 張力なし ---
fprintf('ケースA: 張力なしのシミュレーションを実行中...\n');
params_A = params_base;
params_A.t_start_dynamic = t_start_dynamic;
params_A.t_ramp_dynamic = t_ramp_dynamic;
[t_A, x_A] = ode45(@(t,x) nonlinearForcedODE(t,x,params_A,omega_drive), tspan, x0);
z_A = x_A(:,1);

% --- ケースB: 張力あり ---
fprintf('ケースB: 張力ありのシミュレーションを実行中...\n');
params_B = params_base;
% プリテンションによる剛性増加分を計算
membrane_Area = params_B.w1 * params_B.t1;
strain_y_preload = 2 * z_preload_for_tension^2 / params_B.L^2;
Tension_Force_y = params_B.E1 * strain_y_preload * membrane_Area;
k_tension = 2 * Tension_Force_y / params_B.L;
% 元の剛性にプリテンションによる剛性を加算
params_B.k1_effective = params_B.k1 + k_tension;
params_B.t_start_dynamic = t_start_dynamic;
params_B.t_ramp_dynamic = t_ramp_dynamic;
[t_B, x_B] = ode45(@(t,x) nonlinearForcedODE_tensioned(t,x,params_B,omega_drive), tspan, x0);
z_B = x_B(:,1);

fprintf('---------------------------------------------------\n');
fprintf('比較条件:\n');
fprintf(' プリテンション変位: %.2f mm\n', z_preload_for_tension*1000);
fprintf(' -> これによる剛性増加: %.2f N/m (元のk1=%.2f N/m)\n', k_tension, params_base.k1);
fprintf(' -> 合成剛性 (k1 + k_tension): %.2f N/m\n', params_B.k1_effective);
fprintf(' 両ケースの駆動周波数: %.0f Hz\n', f_drive);
fprintf('---------------------------------------------------\n');

% -------------------------
% 3. 比較プロット
% -------------------------
figure('Name', 'Tension Effect Comparison');
hold on;
plot(t_A, z_A*1000, 'b-', 'DisplayName', 'ケースA: 張力なし');
plot(t_B, z_B*1000, 'r-', 'DisplayName', 'ケースB: 張力あり');
hold off;
title(sprintf('張力の効果の比較 (%.0f Hzで駆動)', f_drive));
xlabel('時間 [s]');
ylabel('変位 z(t) [mm]');
legend;
grid on;
set(gca,'Fontsize',20,'FontWeight','bold');

% -------------------------
% 4. 関数定義
% -------------------------
function p = generate_params(p)
    p.I1 = p.w1 * p.t1^3 / 12; p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.m1 = p.roh1 * p.L * p.w1 * p.t1; p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1);
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3);
    p.m2 = p.S2 * p.w2 * p.roh2; p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2);
end

% 元の（張力なしの）運動方程式
function dx = nonlinearForcedODE(t,x,p,omega)
    z = x(1); dz = x(2);
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2* p.t2) * p.w1;
        F_dynamic = amplitude * cos(omega * time_since_start);
    end
    F_total = F_dynamic;
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end

% 張力あり（剛性が高い）の運動方程式
function dx = nonlinearForcedODE_tensioned(t,x,p,omega)
    z = x(1); dz = x(2);
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2* p.t2) * p.w1;
        F_dynamic = amplitude * cos(omega * time_since_start);
    end
    F_total = F_dynamic;
    denom = p.L^2 - 2*z^2;
    % ★ ばね定数に k1_effective を使用
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1_effective*z - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end