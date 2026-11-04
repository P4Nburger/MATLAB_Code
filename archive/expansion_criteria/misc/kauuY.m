% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc;

% --- ★ 目標とするY方向の静的変位（オフセット）をここで設定 ---
% y変位は常に負になるため、必ず負の値を設定してください
y_static_target = -0.05e-3; % [m] 目標とする静的変位 (例: -0.05 mm)
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
% 2. 荷重計算とシミュレーション
% -------------------------
% パラメータから派生値を計算
params = generate_params(params);

% --- ★ 目標のy変位から、必要なz変位と"Z方向"の静的荷重を計算 ---
if y_static_target > 0
    error('y_static_targetは物理的に負の値である必要があります。');
end
% y = -2z^2/L  =>  z = sqrt(-y*L/2)
z_static_target = sqrt(-y_static_target * params.L / 2);

z_eq = z_static_target;
denom_eq = params.L^2 - 2*z_eq^2;
F_static_required = params.k1*z_eq + (4*params.k2/denom_eq)*z_eq^3;
params.F_static = F_static_required; % "Z方向"にかける静的荷重
params.z_static_target = z_static_target;
% -------------------------------------------------------------

% 初期条件を設定
x0 = [z_static_target; 0]; % 釣り合い点からスタート

fprintf('---------------------------------------------------\n');
fprintf('目標Y変位 %.4f mm を達成するため、\n', y_static_target*1000);
fprintf(' Z変位を %.4f mm に維持します。\n', z_static_target*1000);
fprintf(' そのために必要な【Z方向】の静的荷重: %.4f N\n', F_static_required);
fprintf('---------------------------------------------------\n');

% 周波数スイープを実行
amp = freq_sweep(params, freq_range, tspan, x0);
[~, idx] = max(amp);
f_res = freq_range(idx);
fprintf('シミュレーション結果:\n');
fprintf(' 固有振動数 (オフセット周り): %.2f Hz\n', f_res);

% 共振周波数での時間応答を計算
omega_res = 2*pi*f_res;
[t, x] = ode45(@(t,x) nonlinearForcedODE(t,x,params,omega_res), tspan, x0);
z = x(:,1);
y = -2*z.^2 / params.L;

% -------------------------
% 3. プロット
% -------------------------
figure;
plot(freq_range, amp*1000, 'b', 'LineWidth', 2);
title('振幅応答'); xlabel('周波数 [Hz]'); ylabel('振幅 [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

figure;
plot(t, z*1000, 'b-');
hold on;
yline(z_static_target * 1000, 'r--', 'LineWidth', 1.5, 'Label', 'Z方向の振動中心');
hold off;
title('時間応答 : z方向'); xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;

figure;
plot(t, y*1000, 'b-');
hold on;
yline(y_static_target * 1000, 'r--', 'LineWidth', 1.5, 'Label', 'Y方向の振動中心');
hold off;
title('時間応答 : y方向'); xlabel('時間 [s]'); ylabel('変位 y(t) [mm]');
set(gca,'Fontsize',20,'FontWeight','bold'); grid on;


% -------------------------
% 4. 関数定義 (変更なし)
% -------------------------
function p = generate_params(p)
    p.I1 = p.w1 * p.t1^3 / 12; p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.m1 = p.roh1 * p.L * p.w1 * p.t1; p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1);
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3);
    p.m2 = p.S2 * p.w2 * p.roh2; p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2);
end
function amp = freq_sweep(p, freq_range, tspan, x0)
    amp = zeros(size(freq_range));
    for i = 1:length(freq_range)
        omega = 2*pi*freq_range(i);
        [~,x] = ode45(@(t,x) nonlinearForcedODE(t,x,p,omega), tspan, x0);
        z = x(:,1);
        z_offset = z - p.z_static_target;
        amp(i) = max(abs(z_offset(end-100:end)));
    end
end
function dx = nonlinearForcedODE(t,x,p,omega)
    z = x(1);
    dz = x(2);
    F_dynamic = p.P0 * (p.L - 2* p.t2) * p.w1 * cos(omega*t);
    F_total = F_dynamic + p.F_static; % F_static は Z方向の力
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z ...
           - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end