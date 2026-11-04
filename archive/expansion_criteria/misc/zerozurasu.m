% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ 新しいつり合い点（オフセット）とタイミングを設定 ★ ---
z_offset = 1.0e-3;      % [m] ★ バネの新しいつり合い点の位置 (例: 1.0 mm)
t_start_dynamic = 3.0;    % [s] 音波（動的荷重）をかけ"始める"時間
t_ramp_dynamic = 1.0;     % [s] 動的荷重をこの時間かけて最大振幅にする
% -----------------------------------------------------------

tspan = [0 8];
freq_range = 10:1:200;
params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', 2, 'E1', 197e9, ...
        'w1', 10e-3, 't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.01, ...
        'w2', 7.5e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, ...
        'roh2', 1.24e3, 'k2', 177.035, 'delta2', 0.0278 ...
    );

% -------------------------
% 2. 準備とシミュレーション実行
% -------------------------
params = generate_params(params);

% --- シミュレーションに必要なパラメータを構造体に格納 ---
params.z_offset = z_offset; % ★ 新しいつり合い点をパラメータに追加
params.t_start_dynamic = t_start_dynamic;
params.t_ramp_dynamic = t_ramp_dynamic;
% ---------------------------------------------------

x0 = [0; 0]; % z=0からスタート

fprintf('「つり合い点オフセット」モデルのシミュレーションを開始します...\n');
fprintf('つり合い点: %.2f mm\n', z_offset*1000);
fprintf('まずは、共振周波数を探索します...\n');

% --- ステップ1: 周波数スイープで共振周波数を探す ---
amp = freq_sweep(params, freq_range, tspan, x0);
[~, idx] = max(amp);
f_res = freq_range(idx);
omega_res = 2 * pi * f_res;
fprintf(' -> 共振周波数が見つかりました: %.2f Hz\n', f_res);

% --- ステップ2: 見つかった共振周波数で本番のシミュレーションを実行 ---
[t, x] = ode45(@(t,x) nonlinearForcedODE_offset(t,x,params,omega_res), tspan, x0);
z = x(:,1);
y = -2*z.^2 / params.L;

% -------------------------
% 3. プロット
% -------------------------
figure('Name', 'Frequency Response');
plot(freq_range, amp*1000, 'b-', 'LineWidth', 2);
title('振幅応答'); xlabel('周波数 [Hz]'); ylabel('振幅 [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

figure('Name', 'Z-Direction Displacement');
plot(t, z*1000, 'b-');
hold on;
yline(z_offset * 1000, 'r--', 'LineWidth', 1.5, 'Label', '新しいつり合い点');
xline(t_start_dynamic, 'm--', 'LineWidth', 1.5, 'Label', '音波 開始');
hold off;
title('時間応答 : z方向'); xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

figure('Name', 'Y-Direction Displacement');
plot(t, y*1000, 'b-');
title('時間応答 : y方向'); xlabel('時間 [s]'); ylabel('変位 y(t) [mm]');
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

function amp = freq_sweep(p, freq_range, tspan, x0)
    amp = zeros(size(freq_range));
    h_waitbar = waitbar(0, '周波数スイープを実行中...');
    for i = 1:length(freq_range)
        waitbar(i/length(freq_range), h_waitbar, sprintf('周波数: %d Hz', freq_range(i)));
        omega = 2*pi*freq_range(i);
        [t_sweep,x_sweep] = ode45(@(t,x) nonlinearForcedODE_offset(t,x,p,omega), tspan, x0);
        z_sweep = x_sweep(:,1);
        final_segment = z_sweep(t_sweep > tspan(2)*0.75);
        % 振幅は新しいつり合い点からの差分で評価
        amp(i) = max(abs(final_segment - p.z_offset));
    end
    close(h_waitbar);
end

% ★★★ つり合い点が z_offset になるように改造した運動方程式 ★★★
function dx = nonlinearForcedODE_offset(t,x,p,omega)
    z = x(1);
    dz = x(2);
    
    % ★ 変位を新しいつり合い点からの差分として定義
    z_from_offset = z - p.z_offset;
    
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2* p.t2) * p.w1;
        F_dynamic = amplitude * cos(omega * time_since_start);
    end
    
    F_total = F_dynamic; % 静的な外力はゼロ

    % ★ 復元力の計算に変位の差分 z_from_offset を使用
    denom = p.L^2 - 2*z_from_offset^2;
    ddz = (1 / (p.m1 + 8*p.m2*z_from_offset^2/denom)) * ...
         ( - (p.c1 + 8*p.c2*z_from_offset^2/denom)*dz ...
           - p.k1 * z_from_offset ...
           - (8*p.m2*z_from_offset/denom)*dz^2 ...
           - (4*p.k2/denom)*z_from_offset^3 + F_total );
           
    dx = [dz; ddz];
end