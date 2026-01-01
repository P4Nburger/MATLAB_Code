% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ 比較実験の条件をここで設定 ★ ---
z_preload_for_tension = 0.08e-3; % [m]
f_drive = 120;            % [Hz]
P0_drive = 10;            % [Pa]
t_start_dynamic = 1.0;    % [s]
t_ramp_dynamic = 1.0;     % [s]
% ----------------------------------------------------

tspan = [0 8];
freq_range = 10:1:250;
params = struct( ...
        'L', 30e-3, 'L_side', 60e-3, 'P0', P0_drive, 'E1', 197e9, ...
        'w1', 10e-3, 't1', 0.02e-3, 'roh1', 7.93e3, 'delta1', 0.01, ...
        'w2', 7.5e-3, 'L_up', 15e-3, 't2', 1e-3, 'S2', 126.87e-6, ...
        'roh2', 1.24e3, 'k2', 177.035, 'delta2', 0.0278 ...
    );

% -------------------------
% 2. シミュレーション実行
% -------------------------
params_base = generate_params(params);
x0 = [0;0];

% --- ケースA: プリテンションなし ---
fprintf('ケースA（プリテンションなし）のシミュレーションを開始します...\n');
params_A = params_base;
params_A.t_start_dynamic = t_start_dynamic;
params_A.t_ramp_dynamic = t_ramp_dynamic;
fprintf(' -> ステップ1: 共振周波数を探索中...\n');
amp_A = freq_sweep(params_A, freq_range, tspan, x0, @nonlinearForcedODE_no_tension);
[~, idx_A] = max(amp_A);
f_res_A = freq_range(idx_A);
fprintf(' -> 共振周波数: %.2f Hz\n', f_res_A);
fprintf(' -> ステップ2: 時間応答を計算中...\n');
omega_res_A = 2 * pi * f_res_A;
[t_A, x_A] = ode45(@(t,x) nonlinearForcedODE_no_tension(t,x,params_A,omega_res_A), tspan, x0);
z_A = x_A(:,1);
y_A = -2*z_A.^2 / params_A.L; % ★ Y変位を計算

% --- ケースB: プリテンションあり ---
fprintf('\nケースB（プリテンションあり）のシミュレーションを開始します...\n');
params_B = params_base;
membrane_Area = params_B.w1 * params_B.t1;
strain_y_preload = 2 * z_preload_for_tension^2 / params_B.L^2;
Tension_Force_y = params_B.E1 * strain_y_preload * membrane_Area;
k_tension = 2 * Tension_Force_y / params_B.L;
params_B.k1_effective = params_B.k1 + k_tension;
params_B.t_start_dynamic = t_start_dynamic;
params_B.t_ramp_dynamic = t_ramp_dynamic;
fprintf(' -> ステップ1: 共振周波数を探索中...\n');
amp_B = freq_sweep(params_B, freq_range, tspan, x0, @nonlinearForcedODE_with_tension);
[~, idx_B] = max(amp_B);
f_res_B = freq_range(idx_B);
fprintf(' -> 共振周波数: %.2f Hz\n', f_res_B);
fprintf(' -> ステップ2: 時間応答を計算中...\n');
omega_res_B = 2 * pi * f_res_B;
[t_B, x_B] = ode45(@(t,x) nonlinearForcedODE_with_tension(t,x,params_B,omega_res_B), tspan, x0);
z_B = x_B(:,1);
y_B = -2*z_B.^2 / params_B.L; % ★ Y変位を計算

fprintf('---------------------------------------------------\n');
fprintf('比較結果:\n');
fprintf(' ケースA (プリテンションなし) の剛性 k1: %.2f N/m, 共振周波数: %.2f Hz\n', params_A.k1, f_res_A);
fprintf(' ケースB (プリテンションあり) の剛性 k1_eff: %.2f N/m, 共振周波数: %.2f Hz\n', params_B.k1_effective, f_res_B);
fprintf('---------------------------------------------------\n');

% -------------------------
% 3. 比較プロット
% -------------------------
figure('Name', 'Frequency Response Comparison');
hold on;
plot(freq_range, amp_A*1000, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('ケースA: プリテンションなし (共振点 %.2f Hz)', f_res_A));
plot(freq_range, amp_B*1000, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('ケースB: プリテンションあり (共振点 %.2f Hz)', f_res_B));
plot(f_res_A, max(amp_A)*1000, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
plot(f_res_B, max(amp_B)*1000, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
title('周波数応答の比較'); xlabel('周波数 [Hz]'); ylabel('振幅 [mm]');
legend('Location', 'northwest'); grid on; set(gca,'Fontsize',16,'FontWeight','bold');

figure('Name', 'Z-Direction Time Response Comparison');
hold on;
plot(t_A, z_A*1000, 'b-', 'DisplayName', sprintf('ケースA (%.2f Hzで駆動)', f_res_A));
plot(t_B, z_B*1000, 'r-', 'DisplayName', sprintf('ケースB (%.2f Hzで駆動)', f_res_B));
hold off;
title('時間応答の比較 : z方向');
xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
legend; grid on; set(gca,'Fontsize',16,'FontWeight','bold');

% ★★★ Y方向の時間応答比較プロットを追加 ★★★
figure('Name', 'Y-Direction Time Response Comparison');
hold on;
plot(t_A, y_A*1000, 'b-', 'DisplayName', 'ケースA: プリテンションなし');
plot(t_B, y_B*1000, 'r-', 'DisplayName', 'ケースB: プリテンションあり');
hold off;
title('時間応答の比較 : y方向');
xlabel('時間 [s]');
ylabel('変位 y(t) [mm]');
legend;
grid on;
set(gca,'Fontsize',16,'FontWeight','bold');
% ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★

% -------------------------
% 4. 関数定義 (変更なし)
% -------------------------
function p = generate_params(p)
    p.I1 = p.w1 * p.t1^3 / 12; p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.m1 = p.roh1 * p.L * p.w1 * p.t1; p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1);
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3);
    p.m2 = p.S2 * p.w2 * p.roh2; p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2);
end

function amp = freq_sweep(p, freq_range, tspan, x0, ode_func_handle)
    amp = zeros(size(freq_range));
    h_waitbar = waitbar(0, '周波数スイープを実行中...');
    for i = 1:length(freq_range)
        waitbar(i/length(freq_range), h_waitbar, sprintf('周波数: %d Hz', freq_range(i)));
        omega = 2*pi*freq_range(i);
        [~,x] = ode45(@(t,x) ode_func_handle(t,x,p,omega), tspan, x0);
        z = x(:,1);
        amp(i) = max(abs(z(end-100:end)));
    end
    close(h_waitbar);
end

function dx = nonlinearForcedODE_no_tension(t,x,p,omega)
    z = x(1); dz = x(2);
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2* p.t2) * p.w1;
        F_dynamic = amplitude * cos(omega * time_since_start);
    end
    F_total = F_dynamic;
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end

function dx = nonlinearForcedODE_with_tension(t,x,p,omega)
    z = x(1); dz = x(2);
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2* p.t2) * p.w1;
        F_dynamic = amplitude * cos(omega * time_since_start);
    end
    F_total = F_dynamic;
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1_effective*z - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end