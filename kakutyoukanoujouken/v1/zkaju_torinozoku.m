% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ 静的荷重のタイミングと目標変位を設定 ---
t_ramp_static = 1.0;      % [s] ★ この時間をかけて静的荷重を最大値まで増加させる
z_static_target = 1.0e-3; % [m] ★ 目標とするZ方向の静的変位
% ----------------------------------------------------

tspan = [0 8];
params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', 0, 'E1', 197e9, ...
        'w1', 10e-3, 't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.01, ...
        'w2', 7.5e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, ...
        'roh2', 1.24e3, 'k2', 177.035, 'delta2', 0.0278 ...
    );

% -------------------------
% 2. 荷重計算とシミュレーション
% -------------------------
params = generate_params(params);

z_eq = z_static_target;
denom_eq = params.L^2 - 2*z_eq^2;
F_static_required = params.k1*z_eq + (4*params.k2/denom_eq)*z_eq^3;
params.F_static_max = F_static_required;
params.t_ramp_static = t_ramp_static; % ★ 荷重を取り除く時刻として使用

x0 = [0; 0]; % z=0からスタート

fprintf('シミュレーションを実行します...\n');
fprintf(' %.1f秒かけて荷重を最大にし、その後荷重を取り除きます。\n', t_ramp_static);

[t, x] = ode45(@(t,x) nonlinearForcedODE_release(t,x,params,0), tspan, x0);
z = x(:,1);
y = -2*z.^2 / params.L;

% -------------------------
% 3. プロット
% -------------------------
figure;
plot(t, z*1000, 'b-'); hold on;
yline(z_static_target * 1000, 'k:', 'LineWidth', 1.5, 'Label', '目標変位');
xline(t_ramp_static, 'r--', 'LineWidth', 1.5, 'Label', '荷重除去');
hold off;
title('時間応答 : z方向'); xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

figure;
plot(t, y*1000, 'b-');
title('時間応答 : y方向'); xlabel('時間 [s]'); ylabel('変位 y(t) [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

% --- 荷重履歴のグラフを修正 ---
F_static_history = zeros(size(t));
indices_before_release = t <= params.t_ramp_static;
F_static_history(indices_before_release) = params.F_static_max * (t(indices_before_release) / params.t_ramp_static);
figure;
plot(t, F_static_history, 'm-', 'LineWidth', 2.5);
title('印加される荷重の時間変化'); xlabel('時間 [s]'); ylabel('荷重 [N]');
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

% ★★★ 荷重の印加と除去を行う運動方程式 ★★★
function dx = nonlinearForcedODE_release(t,x,p,omega)
    z = x(1);
    dz = x(2);
    
    F_total = 0;

    if t < p.t_ramp_static
        % --- フェーズ1: 静的荷重をかけて膜をたわませる ---
        ramp_static_factor = t / p.t_ramp_static;
        F_total = p.F_static_max * ramp_static_factor;
    else
        % --- フェーズ2: 静的荷重をゼロにする ---
        F_total = 0; 
    end
    
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z ...
           - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
           
    dx = [dz; ddz];
end