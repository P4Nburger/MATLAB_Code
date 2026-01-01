% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc;

% --- ★ ご希望の現象が見やすいように設定 ---
z0_initial   = 2e-3;    % 初期変位を大きく設定 (2mm)
dz0_initial  = 0;       % 初期速度はゼロ
tspan = [0 15];         % ★ シミュレーション時間を長くする (15秒)
% -----------------------------------------

freq_range = 10:1:200;

params = struct( ...
        'L', 30e-3, ...
        'P0', 2, ...
        'E1', 197e9, ...
        'w1', 10e-3, ...
        't1', 0.02e-3, ...
        'roh1', 7.93e3, ...
        'delta1', 0.004, ... % ★ 減衰を意図的に小さくして、移行をゆっくりにする
        'w2', 7.5e-3, ...
        'k2', 177.035, ...
        'delta2', 0.0278, ...
        'L_side', 60e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, 'roh2', 1.24e3 ...
    );

x0 = [z0_initial; dz0_initial];

% -------------------------
% 3. シミュレーション実行とプロット
% -------------------------
params = generate_params(params);
amp = freq_sweep(params, freq_range, tspan, x0);
[~, idx] = max(amp);
f_res = freq_range(idx);
fprintf('共振周波数: %.2f Hz\n', f_res);

% 共振周波数で時間応答を計算
omega_res = 2*pi*f_res;
[t, x] = ode45(@(t,x) nonlinearForcedODE(t,x,params,omega_res), tspan, x0);
z = x(:,1);

% グラフをプロット
figure;
plot(t, z*1000, 'b-', 'LineWidth', 1); % 変位をmm単位でプロット
hold on;

% --- 補助線を追加して現象を分かりやすくする ---
% 初期変位のライン
yline(z0_initial*1000, 'r--', 'LineWidth', 1.5, 'Label', '初期変位');
% 最終的な振動中心のライン
yline(0, 'k--', 'LineWidth', 1, 'Label', '最終的な振動中心 (z=0)');
% 振動中心の移動の目安となる包絡線 (理論的な減衰曲線)
zeta = params.delta1 / (2*pi); % 減衰比 (概算)
omega_n = sqrt(params.k1 / params.m1); % 固有角周波数 (概算)
envelope = z0_initial * exp(-zeta * omega_n * t);
plot(t, envelope*1000, 'g:', 'LineWidth', 2);
plot(t, -envelope*1000, 'g:', 'LineWidth', 2);
hold off;

set(gca,'Fontsize',16,'FontWeight','bold');
xlabel('時間 [s]');
ylabel('変位 z(t) [mm]');
title('振動中心が初期位置からゼロへ移行する様子');
legend('z(t)の変位', '初期変位', '最終的な振動中心', '振動中心の減衰の目安', 'Location', 'northeast');
grid on;


% -----------------------------------------------------
% --- ここから下は関数の定義（必ずファイルの最後に置く）---
% -----------------------------------------------------

function p = generate_params(p)
    p.I1 = p.w1 * p.t1^3 / 12;
    p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.m1 = p.roh1 * p.L * p.w1 * p.t1;
    p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1);
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3);
    p.m2 = p.S2 * p.w2 * p.roh2;
    p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2);
end

function dx = nonlinearForcedODE(t,x,p,omega)
    z = x(1);
    dz = x(2);
    denom = p.L^2 - 2*z^2;
    F = p.P0 * (p.L - 2* p.t2) * p.w1 * cos(omega*t);
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z ...
           - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F );
    dx = [dz; ddz];
end

function amp = freq_sweep(p, freq_range, tspan, x0)
    amp = zeros(size(freq_range));
    for i = 1:length(freq_range)
        omega = 2*pi*freq_range(i);
        [~,x] = ode45(@(t,x) nonlinearForcedODE(t,x,p,omega), tspan, x0);
        z = x(:,1);
        amp(i) = max(abs(z(round(end/2):end))); % 後半の振幅を取得
    end
end