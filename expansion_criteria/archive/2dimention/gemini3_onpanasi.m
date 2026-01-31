% =========================================================================
% 2自由度 (2-DOF) 拡張メタマテリアル シミュレーション
% 【音波なし】: 初期たわみによる「フレームの収縮（マイナス変位）」を確認
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
tspan = [0 1];       % [s] シミュレーション時間

% --- 目標とする状態の設定 ---
z_target_initial = 1.0e-3; % [m] 初期たわみ量 (1mm)

% --- 音波の設定 (無効化) ---
f_drive = 120;       
P0_drive = 0;        % [Pa] 音圧0

% --- 物理定数 ---
params = struct();
params.L = 30e-3;      
params.w1 = 10e-3;     
params.t1 = 0.02e-3;   
params.E1 = 197e9;     
params.rho1 = 7930;    
params.S2 = 126.87e-6; 
params.rho2 = 1240;    
params.k2 = 177.0;     
params.delta = 0.05;   
params.P0 = P0_drive;
params.freq = f_drive;

% -------------------------------------------------------------------------
% 2. 派生パラメータと平衡点計算
% -------------------------------------------------------------------------
params.m1 = params.rho1 * params.L * params.w1 * params.t1;
params.m2 = params.rho2 * params.S2 * params.L;
params.I1 = params.w1 * params.t1^3 / 12;
params.k1 = 384 * params.E1 * params.I1 / (params.L^3);
params.K_couple = params.E1 * (params.w1 * params.t1) / params.L;
params.c1 = 2 * params.delta * sqrt(params.m1 * params.k1);
params.c2 = 2 * params.delta * sqrt(params.m2 * params.k2);

% --- 平衡点の計算 ---
z = z_target_initial;
term = (params.k1 * z * params.L) / (4 * z * params.K_couple);
y_eq = term + (2 * z^2 / params.L); % これは「収縮量（大きさ）」なのでプラスの値

term_couple = params.K_couple * (y_eq - 2 * z^2 / params.L);
y_natural = y_eq + term_couple / params.k2;
params.y_natural = y_natural;

fprintf('--- 設定完了 ---\n');
fprintf('目標初期たわみ (Z): %.3f mm\n', z*1000);
fprintf('平衡点 収縮量 (Y): %.3f mm\n', y_eq*1000);

% -------------------------------------------------------------------------
% 3. シミュレーション実行
% -------------------------------------------------------------------------
x0 = [z_target_initial; y_eq; 0; 0]; 

fprintf('シミュレーション実行中...\n');
options = odeset('RelTol',1e-6, 'AbsTol',1e-9);
[t, x] = ode45(@(t,x) equations_2DOF_NoSound(t, x, params), tspan, x0, options);

z_res = x(:,1);
y_res = x(:,2);

% -------------------------------------------------------------------------
% 4. プロット (Y方向をマイナス表示に変更)
% -------------------------------------------------------------------------
figure('Position', [100, 100, 800, 800], 'Color', 'w');

% (1) Z方向変位 (膜)
subplot(2,1,1);
plot(t, z_res*1000, 'b-', 'LineWidth', 1.5);
hold on;
yline(z_target_initial*1000, 'g--', 'LineWidth', 1.5, 'Label', '初期たわみ位置');
yline(0, 'k-', 'LineWidth', 1); % 0ライン
hold off;
ylabel('Z方向変位 [mm]', 'FontSize', 12);
title('膜の変位 (Z) - 初期たわみ有り', 'FontSize', 14);
grid on;
ylim([0, max(z_res)*1000*1.2]);

% (2) Y方向変位 (フレーム) - マイナス方向＝押しつぶされた状態
subplot(2,1,2);

% ★★★ 修正箇所：符号を反転させてプロット ★★★
% y_res は「収縮量（大きさ）」なので、座標としては「マイナス」にする
y_plot = -y_res * 1000; 
y_target_plot = -y_eq * 1000;

plot(t, y_plot, 'r-', 'LineWidth', 1.5);
hold on;
yline(0, 'k-', 'LineWidth', 1, 'Label', '変形なし(0)'); % 基準線
yline(y_target_plot, 'g--', 'LineWidth', 1.5, 'Label', '初期たわみによる収縮位置');
hold off;

ylabel('Y方向変位 [mm]', 'FontSize', 12);
title('フレームの変位 (Y) - 負の方向に押しつぶされている', 'FontSize', 14);
grid on;
xlabel('時間 [s]', 'FontSize', 12);

% グラフの範囲調整（0とマイナス値が入るように）
ylim([y_target_plot*1.2, 0.01]); 


% =========================================================================
% 2自由度 運動方程式 (音波なし)
% =========================================================================
function dx = equations_2DOF_NoSound(t, x, p)
    z  = x(1);
    y  = x(2);
    dz = x(3);
    dy = x(4);
    
    coupling_term = y - (2 * z^2 / p.L);
    
    F_restore_z = p.k1 * z - p.K_couple * coupling_term * (4 * z / p.L);
    F_restore_y = p.k2 * (y - p.y_natural) + p.K_couple * coupling_term;
    
    F_sound = 0;
    
    F_damping_z = p.c1 * dz;
    F_damping_y = p.c2 * dy;
    
    ddz = (F_sound - F_damping_z - F_restore_z) / p.m1;
    ddy = (0 - F_damping_y - F_restore_y) / p.m2;
    
    dx = [dz; dy; ddz; ddy];
end