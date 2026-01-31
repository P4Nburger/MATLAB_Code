% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ シミュレーションの動作をここで設定 ★ ---
z_static_target = 1.0e-3; % [m] 目標とする静的変位
t_ramp_static = 1.0;      % [s] この時間で静的変位に到達させる
% -----------------------------------------------------------

tspan = [0 8]; % 動きを確認するためシミュレーション時間を短めに設定
params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', 0, 'E1', 197e9, ...
        'w1', 10e-3, 't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.01, ...
        'w2', 7.5e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, ...
        'roh2', 1.24e3, 'k2', 177.035, 'delta2', 0.0278 ...
    );

% -------------------------
% 2. 準備とシミュレーション実行
% -------------------------
params = generate_params(params);

% --- ステップ1: 静的荷重の計算 ---
z_eq = z_static_target;
denom_eq = params.L^2 - 2*z_eq^2;
F_static_required = params.k1*z_eq + (4*params.k2/denom_eq)*z_eq^3;
params.F_static_max = F_static_required;
params.t_ramp_static = t_ramp_static;

% --- ステップ2: プリテンションによる剛性増加分を計算 ---
membrane_Area = params.w1 * params.t1;
strain_y_preload = 2 * z_static_target^2 / params.L^2;
Tension_Force_y = params.E1 * strain_y_preload * membrane_Area;
k_tension = 2 * Tension_Force_y / params.L;
params.k1_effective = params.k1 + k_tension;

x0 = [0; 0];

fprintf('新しいプリテンションモデルのシミュレーションを開始します...\n');

[t, x] = ode45(@(t,x) nonlinearForcedODE_staged(t,x,params,0), tspan, x0);
z = x(:,1);
y = -2*z.^2 / params.L;

% -------------------------
% 3. プロット
% -------------------------
figure('Name', 'Z-Direction Displacement');
plot(t, z*1000, 'b-', 'LineWidth', 2);
hold on;
xline(t_ramp_static, 'r--', 'LineWidth', 1.5, 'Label', 'プリテンション ON');
yline(z_static_target * 1000, 'k:', 'Label', '目標変位');
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

% ★★★ 時間経過でモデルが切り替わる運動方程式 ★★★
function dx = nonlinearForcedODE_staged(t,x,p,omega)
    z = x(1);
    dz = x(2);
    
    F_total = 0;
    k_current = p.k1; % デフォルトは元の剛性

    if t < p.t_ramp_static
        % --- フェーズ1: 静的荷重をランプアップ ---
        ramp_static_factor = t / p.t_ramp_static;
        F_total = p.F_static_max * ramp_static_factor;
        k_current = p.k1; % この間はまだ元の剛性
        
    else
        % --- フェーズ2: プリテンション効果をON ---
        F_total = 0; % 外力としての静的荷重はゼロにする
        k_current = p.k1_effective*100; % 剛性を高める
    end
    
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (p.c1 + 8*p.c2*z^2/denom)*dz - k_current*z ... % ★ 現在の剛性を使用
           - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
           
    dx = [dz; ddz];
end