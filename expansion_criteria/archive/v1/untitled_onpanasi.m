% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ 静的荷重のタイミングと目標変位を設定 ---
t_ramp_static = 1.0;     % [s] この時間をかけて静的荷重を最大値まで増加させる
y_static_target = -0.05e-3; % [m] 目標とするY方向の静的変位
% ----------------------------------------------------

tspan = [0 8];
params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', 0, 'E1', 197e9, ... % ★ P0 (音圧) を 0 に設定
        'w1', 10e-3, 't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.01, ...
        'w2', 7.5e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, ...
        'roh2', 1.24e3, 'k2', 177.035, 'delta2', 0.0278 ...
    );

% -------------------------
% 2. 荷重計算とシミュレーション
% -------------------------
params = generate_params(params);

if y_static_target > 0
    error('y_static_targetは物理的に負の値である必要があります。');
end
z_static_target = sqrt(-y_static_target * params.L / 2);
z_eq = z_static_target;
denom_eq = params.L^2 - 2*z_eq^2;
F_static_required = params.k1*z_eq + (4*params.k2/denom_eq)*z_eq^3;
params.F_static_max = F_static_required;
params.t_ramp_static = t_ramp_static;
% ★ 動的荷重関連のパラメータは不要
% params.t_start_dynamic = t_start_dynamic;
% params.t_ramp_dynamic = t_ramp_dynamic;
params.z_static_target = z_static_target;

x0 = [0; 0]; % z=0からスタート

fprintf('---------------------------------------------------\n');
fprintf('目標Y変位 %.4f mm\n', y_static_target*1000);
fprintf(' そのために必要な最大静的荷重(Z方向): %.4f N\n', F_static_required);
fprintf(' 静的荷重 ramp時間: %.1f 秒\n', t_ramp_static);
fprintf('---------------------------------------------------\n');

% ★ 周波数スイープは不要なので削除

fprintf('シミュレーションを実行します...\n');
% ★ ode45を直接実行 (omegaは使われないのでダミーで0を渡す)
[t, x] = ode45(@(t,x) nonlinearForcedODE(t,x,params,0), tspan, x0);
z = x(:,1);
y = -2*z.^2 / params.L;

% -------------------------
% 3. プロット
% -------------------------
% ★ 周波数応答グラフは不要なので削除

% 時間応答グラフ (z方向)
figure;
plot(t, z*1000, 'b-'); hold on;
yline(z_static_target * 1000, 'r--', 'LineWidth', 1.5, 'Label', '目標変位'); hold off;
title('時間応答 : z方向'); xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

% 時間応答グラフ (y方向)
figure;
plot(t, y*1000, 'b-'); hold on;
yline(y_static_target * 1000, 'r--', 'LineWidth', 1.5, 'Label', '目標変位'); hold off;
title('時間応答 : y方向'); xlabel('時間 [s]'); ylabel('変位 y(t) [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

% 静的荷重の履歴グラフ
F_static_history = params.F_static_max * min(1.0, t / params.t_ramp_static);
figure;
plot(t, F_static_history, 'm-', 'LineWidth', 2.5, 'DisplayName', '静的荷重');
title('印加される荷重の時間変化'); xlabel('時間 [s]'); ylabel('荷重 [N]');
legend;
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

% ★ 周波数スイープ関数は不要なので削除

% ★ 運動方程式を簡略化
function dx = nonlinearForcedODE(t,x,p,omega)
    z = x(1);
    dz = x(2);
    
    % 静的荷重のみを計算
    ramp_static_factor = min(1.0, t / p.t_ramp_static);
    F_static_current = p.F_static_max * ramp_static_factor;
    
    F_total = F_static_current; % 動的荷重はゼロ
    
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z ...
           - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end