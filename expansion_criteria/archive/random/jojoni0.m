% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc;

% --- 基本設定 ---
z0_initial   = 2e-3;    % 初期変位 (2mm)
dz0_initial  = 0;       % 初期速度はゼロ
tspan = [0 20];         % 変化が分かりやすいよう時間を少し延ばす
freq_range = 10:1:200;

% --- ★ 振動に応じて静的力を減衰させるための追加パラメータ ---
% これらの値を調整すると、中心の移動の仕方が変わります
alpha = 0.8;  % 振幅検出の感度 (大きいと速く反応)
gamma = 0.3;  % 振幅検出の時定数 (大きいとAが速く減衰)
beta  = 500;  % 静的力が振幅Aに応じて減衰する速さ (大きいと速く中心が移動)
% -----------------------------------------------------------

params = struct( ...
        'L', 30e-3, 'P0', 2, 'E1', 197e9, 'w1', 10e-3, 't1', 0.02e-3, ...
        'roh1', 7.93e3, 'delta1', 0.02, 'w2', 7.5e-3, 'k2', 177.035, ...
        'delta2', 0.0278, 'L_side', 60e-3, 'L_up', 15e-3, 't2', 1e-3, ...
        'S2', 126.87e-6, 'roh2', 1.24e3, ...
        'alpha', alpha, 'gamma', gamma, 'beta', beta ... % ★ パラメータに追加
    );

% ★ 初期条件ベクトルを3次元に拡張 [z; dz; A]
x0 = [z0_initial; dz0_initial; 0]; % A(0) = 0 からスタート

% -------------------------
% 2. シミュレーション実行とプロット
% -------------------------
params = generate_params(params);

% 初期位置(z0)で釣り合うために必要な静的な初期力を計算
z_init = z0_initial;
denom = params.L^2 - 2*z_init^2;
restoring_force_at_z0 = - (params.k1*z_init + (4*params.k2/denom)*z_init^3);
params.F_static_initial = -restoring_force_at_z0;

% 周波数スイープを実行して共振周波数を探す
amp = freq_sweep(params, freq_range, tspan, x0);
[~, idx] = max(amp);
f_res = freq_range(idx);
fprintf('共振周波数: %.2f Hz\n', f_res);

% 共振周波数で本計算
omega_res = 2*pi*f_res;
[t, x] = ode45(@(t,x) nonlinearForcedODE(t,x,params,omega_res), tspan, x0);
z_sim = x(:,1); % z変位
A_sim = x(:,3); % 振幅インジケータA

% グラフをプロット
figure;
plot(t, z_sim*1000, 'b-');
hold on;

% 振動中心の移動ラインを計算してプロット
center_of_vibration = z0_initial * exp(-params.beta * A_sim);
plot(t, center_of_vibration*1000, 'r--', 'LineWidth', 2.5);
yline(0, 'k--');
hold off;

set(gca,'Fontsize',16,'FontWeight','bold');
xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
title('【最終版】振動の大きさに応じて振動中心が移動する様子');
legend('z(t)の変位', '振動中心の移動', '最終的な振動中心');
grid on;

% -----------------------------------------------------
% --- ここから下は関数の定義 ---
% -----------------------------------------------------

function p = generate_params(p)
    p.I1 = p.w1 * p.t1^3 / 12; p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.m1 = p.roh1 * p.L * p.w1 * p.t1; p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1);
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3);
    p.m2 = p.S2 * p.w2 * p.roh2; p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2);
end

% ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
% 運動方程式を改造: 状態変数を3つ [z, dz, A] にし、Aに応じて静的力が変化
% ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
function dx = nonlinearForcedODE(t,x,p,omega)
    % 状態変数を分離
    z = x(1); % 変位
    dz = x(2); % 速度
    A = x(3); % 振幅インジケータ
    
    % --- 力を計算 ---
    % 1. 振幅インジケータAに応じて減衰する静的な力
    F_static = p.F_static_initial * exp(-p.beta * A);
    
    % 2. 元々の外力（強制振動）
    F_dynamic = p.P0 * (p.L - 2* p.t2) * p.w1 * cos(omega*t);
    
    % 3. 全ての力を合算
    F_total = F_dynamic + F_static;
    
    % --- 微分方程式 ---
    % 1. zの2階微分（加速度）
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z ...
           - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
           
    % 2. Aの1階微分（振幅インジケータの変化率）
    dA = p.alpha * abs(z) - p.gamma * A;
    
    % 3つの状態変数の微分をベクトルとして返す
    dx = [dz; ddz; dA];
end
% ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★

function amp = freq_sweep(p, freq_range, tspan, x0)
    amp = zeros(size(freq_range));
    for i = 1:length(freq_range)
        omega = 2*pi*freq_range(i);
        [~,x] = ode45(@(t,x) nonlinearForcedODE(t,x,p,omega), tspan, x0);
        z = x(:,1);
        amp(i) = max(abs(z(round(end/2):end)));
    end
end