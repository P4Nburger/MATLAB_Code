% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ シミュレーションの動作をここで設定 ---
z_static_target = 4.8e-3; % [m] まず移動させる安定点 (約4.8mm)
t_ramp_static = 1.0;      % [s] 安定点に到達させる時間
t_start_dynamic = 3.0;    % [s] ★「飛び移り」を誘発する音波をかけ始める時間
P0_kick = 140;             % [Pa] ★ 飛び移りを起こすための、大きな音圧振幅
f_dynamic = 2;            % [Hz] ★ 低い周波数でゆっくりと大きな力をかける
% ----------------------------------------------------

tspan = [0 8]; % シミュレーション時間
params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', P0_kick, 'E1', 197e9, ... % P0をP0_kickに変更
        'w1', 10e-3, 't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.1, ... % 減衰を少し上げて見やすくする
        'w2', 7.5e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, ...
        'roh2', 1.24e3, 'k2', 177.035, 'delta2', 0.0278 ...
    );

% -------------------------
% 2. 準備とシミュレーション実行
% -------------------------
params = generate_params(params);

z_eq = z_static_target;
denom_eq = params.L^2 - 2*z_eq^2;
F_static_required = params.k1*z_eq - (4*params.k2/denom_eq)*z_eq^3;
params.F_static_max = F_static_required;

params.t_ramp_static = t_ramp_static;
params.t_start_dynamic = t_start_dynamic;
params.t_ramp_dynamic = 0.1; % 急激に力を加えるためramp時間は短く

x0 = [0; 0];
omega_dynamic = 2 * pi * f_dynamic;

fprintf('「飛び移り」現象のシミュレーションを開始します...\n');

[t, x] = ode45(@(t,x) nonlinearForcedODE_bistable(t,x,params,omega_dynamic), tspan, x0);
z = x(:,1);
y = -2*z.^2 / params.L;

% -------------------------
% 3. プロット
% -------------------------
figure('Name', 'Z-Direction Displacement Snap-Through');
plot(t, z*1000, 'b-', 'LineWidth', 2);
hold on;
xline(t_start_dynamic, 'm--', 'LineWidth', 1.5, 'Label', '大振幅の動的荷重 開始');
yline(0, 'k:'); % z=0 のライン
hold off;
title('時間応答 : z方向（飛び移り現象）');
xlabel('時間 [s]');
ylabel('変位 z(t) [mm]');
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

function dx = nonlinearForcedODE_bistable(t,x,p,omega)
    z = x(1);
    dz = x(2);
    
    ramp_static_factor = min(1.0, t / p.t_ramp_static);
    F_static_current = p.F_static_max * ramp_static_factor;
    
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0; % P0に面積などを掛ける前の値
        % 実際の力に変換
        Force_amplitude = amplitude * (p.L - 2* p.t2) * p.w1;
        F_dynamic = Force_amplitude * cos(omega * time_since_start);
    end

    F_total = F_dynamic + F_static_current;
    
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (p.c1 + 8*p.c2*z^2/denom)*dz ...
           + p.k1*z ...
           - (8*p.m2*z/denom)*dz^2 ...
           - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end