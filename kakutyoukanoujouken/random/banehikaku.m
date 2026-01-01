% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ モデルの種類をここで選択 ★ ---
% 1: 通常のバネ (線形)
% 2: 硬化バネ (z=0から離れるほど硬くなる)
% 3: 軟化バネ (z=0付近が最も硬い)
model_type = 3;
% ------------------------------------

% --- シミュレーション設定 ---
tspan = [0 8];
f_dynamic = 80;
P0_drive = 10;

params = struct( ...
        'L', 30e-3, 'P0', P0_drive, 'E1', 197e9, 'w1', 10e-3, ...
        't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.01 ...
    );
params.m1 = params.roh1 * params.L * params.w1 * params.t1;
params.I1 = params.w1 * params.t1^3 / 12;
params.k1 = 384 * params.E1 * params.I1 / (params.L^3);
params.c1 = params.delta1 / pi * sqrt(params.m1 * params.k1);
params.k3 = params.k1 * 1e6; % 非線形項の強さを決める係数

% --- シミュレーション実行 ---
x0 = [0; 0];
omega_dynamic = 2 * pi * f_dynamic;
[t, x] = ode45(@(t,x) simpleODE(t,x,params,omega_dynamic,model_type), tspan, x0);
z = x(:,1);

% -------------------------
% 2. プロット
% -------------------------
% --- 時間応答グラフ ---
figure('Name', 'Time Response');
plot(t, z*1000, 'b-');
title(sprintf('時間応答 (モデルタイプ %d)', model_type));
xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
grid on; set(gca,'Fontsize',16);

% --- バネの特性グラフ (力 vs 変位) ---
figure('Name', 'Force vs Displacement');
z_range = linspace(-max(abs(z)), max(abs(z)), 200);
F_restore = zeros(size(z_range));
for i = 1:length(z_range)
    F_restore(i) = get_restoring_force(z_range(i), params, model_type);
end
plot(z_range*1000, F_restore, 'r-', 'LineWidth', 3);
title(sprintf('バネの特性（復元力 vs 変位）- モデル %d', model_type));
xlabel('変位 z [mm]'); ylabel('復元力 F [N]');
grid on; set(gca,'Fontsize',16);

% -------------------------
% 3. 関数定義
% -------------------------
% 復元力を計算する部分を独立した関数に
function F_r = get_restoring_force(z, p, type)
    if type == 1 % 通常バネ
        F_r = p.k1 * z;
    elseif type == 2 % 硬化バネ
        F_r = p.k1 * z + p.k3 * z.^3;
    elseif type == 3 % 軟化バネ
        F_r = p.k1 * z - p.k3 * z.^3;
    else
        F_r = p.k1 * z;
    end
end

% ODEソルバー用の関数
function dx = simpleODE(t,x,p,omega,type)
    z = x(1);
    dz = x(2);
    
    % 復元力の計算
    F_restore = get_restoring_force(z, p, type);
    
    % 動的荷重
    F_dynamic = p.P0 * (p.L) * p.w1 * cos(omega*t);
    
    % 運動方程式 (簡略版)
    ddz = (1 / p.m1) * ( -p.c1*dz - F_restore + F_dynamic );
    
    dx = [dz; ddz];
end