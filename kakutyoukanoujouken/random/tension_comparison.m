% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ 比較実験用の設定 ★ ---
f_drive = 120;            % [Hz] ★ 駆動周波数を固定 (例: 120 Hz)
z_static_target_tension = 1.0e-3; % [m] ★「張力あり」ケースの静的変位
% ----------------------------------------------------

t_ramp_static = 1.0;
t_start_dynamic = 1.0; % すぐに振動を開始
t_ramp_dynamic = 1.0;
tspan = [0 8];

params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', 2, 'E1', 197e9, ...
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
fprintf('ケースA: 張力なしでシミュレーション中...\n');
params_A = params_base;
params_A.F_static_max = 0; % 静的荷重なし
params_A.t_ramp_static = t_ramp_static;
params_A.t_start_dynamic = t_start_dynamic;
params_A.t_ramp_dynamic = t_ramp_dynamic;
[t_A, x_A] = ode45(@(t,x) nonlinearForcedODE_timed(t,x,params_A,omega_drive), tspan, x0);
z_A = x_A(:,1);

% --- ケースB: 張力あり ---
fprintf('ケースB: 張力ありでシミュレーション中...\n');
params_B = params_base;
z_eq_B = z_static_target_tension;
denom_eq_B = params_B.L^2 - 2*z_eq_B^2;
F_static_required_B = params_B.k1*z_eq_B + (4*params_B.k2/denom_eq_B)*z_eq_B^3;
params_B.F_static_max = F_static_required_B;
params_B.t_ramp_static = t_ramp_static;
params_B.t_start_dynamic = t_start_dynamic;
params_B.t_ramp_dynamic = t_ramp_dynamic;
[t_B, x_B] = ode45(@(t,x) nonlinearForcedODE_timed(t,x,params_B,omega_drive), tspan, x0);
z_B = x_B(:,1);

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

function dx = nonlinearForcedODE_timed(t,x,p,omega)
    z = x(1); dz = x(2);
    ramp_static_factor = min(1.0, t / p.t_ramp_static);
    F_static_current = p.F_static_max * ramp_static_factor;
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2* p.t2) * p.w1;
        F_dynamic = amplitude * cos(omega * time_since_start);
    end
    F_total = F_dynamic + F_static_current;
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end