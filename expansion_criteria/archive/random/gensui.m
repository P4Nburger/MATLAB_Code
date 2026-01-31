% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★「中心付近の抵抗」の強さと範囲を設定 ★ ---
extra_damping_factor = 80;  % z=0での追加の抵抗の強さ (0で抵抗なし)
damping_width = 0.2e-3;   % [m] 抵抗が強くなるz=0周りの範囲の広さ (例: 0.2 mm)
% ----------------------------------------------------

% --- 動的荷重の設定 ---
t_start_dynamic = 1.0;    % [s] 音波をかけ"始める"時間
t_ramp_dynamic = 2.0;     % [s] 動的荷重をこの時間かけて最大振幅にする
f_dynamic = 80;           % [Hz] かける音波の周波数
% ----------------------------------------------------

tspan = [0 8]; % シミュレーション時間
params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', 10, 'E1', 197e9, ... % P0を少し上げて動きを大きくする
        'w1', 10e-3, 't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.01, ...
        'w2', 7.5e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, ...
        'roh2', 1.24e3, 'k2', 177.035, 'delta2', 0.0278 ...
    );

% -------------------------
% 2. 準備とシミュレーション実行
% -------------------------
params = generate_params(params);

params.t_start_dynamic = t_start_dynamic;
params.t_ramp_dynamic = t_ramp_dynamic;
params.extra_damping_factor = extra_damping_factor;
params.damping_width = damping_width;

x0 = [0; 0];
omega_dynamic = 2 * pi * f_dynamic;

fprintf('「中心付近の抵抗増加」モデルのシミュレーションを開始します...\n');

[t, x] = ode45(@(t,x) nonlinearForcedODE_center_damping(t,x,params,omega_dynamic), tspan, x0);
z = x(:,1);
dz = x(:,2); % ★ 速度データも取得

% -------------------------
% 3. プロット
% -------------------------
figure('Name', 'Z-Direction Displacement');
plot(t, z*1000, 'b-');
hold on;
yline(0, 'k:');
patch([t fliplr(t)], [ones(size(t))*params.damping_width*1000 fliplr(-ones(size(t))*params.damping_width*1000)], ...
      'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '抵抗増加エリア');
hold off;
title('時間応答 : z方向（中心抵抗モデル）');
xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
legend;
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

% ★★★ 現象を可視化するための「位相図」プロット ★★★
figure('Name', 'Phase Portrait');
plot(z*1000, dz*1000, 'b-');
xlabel('変位 z [mm]');
ylabel('速度 dz/dt [mm/s]');
title('位相図（変位 vs 速度）');
grid on;
set(gca,'Fontsize',20,'FontWeight','bold');
hold on;
% 抵抗が増加するエリアを可視化
xline(params.damping_width*1000, 'r--', 'LineWidth', 2);
xline(-params.damping_width*1000, 'r--', 'LineWidth', 2);
hold off;


% -------------------------
% 4. 関数定義
% -------------------------
function p = generate_params(p)
    p.I1 = p.w1 * p.t1^3 / 12; p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.m1 = p.roh1 * p.L * p.w1 * p.t1; p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1);
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3);
    p.m2 = p.S2 * p.w2 * p.roh2; p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2);
end

% ★★★ 中心付近で抵抗（減衰）が大きくなる運動方程式 ★★★
function dx = nonlinearForcedODE_center_damping(t,x,p,omega)
    z = x(1);
    dz = x(2);
    
    % --- 動的荷重 ---
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2* p.t2) * p.w1;
        F_dynamic = amplitude * cos(omega * time_since_start);
    end
    F_total = F_dynamic; % 静的荷重はなし

    % ★★★ z=0 付近で増加する抵抗（減衰）を計算 ★★★
    c1_extra = p.c1 * p.extra_damping_factor * exp(-(z / p.damping_width)^2);
    c1_effective = p.c1 + c1_extra;
    % ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (c1_effective + 8*p.c2*z^2/denom)*dz ... % ★ 実効減衰係数を使用
           - p.k1*z ...
           - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
           
    dx = [dz; ddz];
end