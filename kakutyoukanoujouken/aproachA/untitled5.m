% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ シミュレーションの動作をここで設定 ---
z_static_target = 1.0e-3; % [m] 目標とする「初期たわみ」
t_ramp_static = 1.0;      % [s] この時間で初期たわみを作成
t_start_dynamic = 3.0;    % [s] 音波（動的荷重）をかけ"始める"時間
t_ramp_dynamic = 1.0;     % [s] 動的荷重をこの時間かけて最大振幅にする
f_dynamic = 85;           % [Hz] かける音波の周波数
% ----------------------------------------------------

tspan = [0 8];
params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', 4, 'E1', 197e9, ...
        'w1', 10e-3, 't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.01, ...
        'w2', 7.5e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, ...
        'roh2', 1.24e3, 'k2', 177.035, 'delta2', 0.0278 ...
    );

% -------------------------
% 2. 準備とシミュレーション実行
% -------------------------
params = generate_params(params);

z_eq = z_static_target;
denom_eq = params.L^2 - 2*z_eq^2;
F_static_required = params.k1*z_eq + (4*params.k2/denom_eq)*z_eq^3;
params.F_static_max = F_static_required;

params.t_ramp_static = t_ramp_static;
params.t_start_dynamic = t_start_dynamic;
params.t_ramp_dynamic = t_ramp_dynamic;

x0 = [0; 0];
omega_dynamic = 2 * pi * f_dynamic;

fprintf('「アプローチA：静的荷重モデル」のシミュレーションを開始します...\n');
[t, x] = ode45(@(t,x) nonlinearForcedODE_timed(t,x,params,omega_dynamic), tspan, x0);
z = x(:,1);

% -------------------------
% 3. プロット
% -------------------------
figure('Name', 'Z-Direction Displacement');
plot(t, z*1000, 'b-');
hold on;
yline(z_static_target * 1000, 'r--', 'LineWidth', 1.5, 'Label', '目標のたわみ');
xline(t_start_dynamic, 'm--', 'LineWidth', 1.5, 'Label', '音波 開始');
hold off;
title('時間応答 : z方向'); xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

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
    z = x(1);
    dz = x(2);
    
    % 静的荷重を t_ramp_static かけて徐々に増加させ、その後は最大値で「維持」
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
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z ...
           - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end