% -------------------------
% 1. パラメータ設定
% -------------------------
clear; clc; close all;

% --- ★ シミュレーションの動作をここで設定 ★ ---
z_static_target = 0.5e-3; % [m] ★ 目標とする初期変位（静的たわみ）
t_ramp_static = 1.0;      % [s] ★ 上記の変位にこの時間をかけて到達させる

% ケースBで適用するプリテンションの強さ
z_preload_for_tension = 0.05e-3; % [m] 剛性計算に使う仮想的なZ変位

% 両ケース共通でかける音波の条件
P0_drive = 10;
t_start_dynamic = 3.0;    % [s] 音波をかけ"始める"時間
t_ramp_dynamic = 1.0;     % [s] 動的荷重をこの時間かけて最大振幅にする
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

% ★ 両ケース共通でかける静的荷重を計算
z_eq_target = z_static_target;
denom_eq_target = params_base.L^2 - 2*z_eq_target^2;
F_static_required = params_base.k1*z_eq_target + (4*params_base.k2/denom_eq_target)*z_eq_target^3;

% --- ケースA: プリテンションなし + 静的荷重あり ---
fprintf('ケースA（プリテンションなし）のシミュレーションを開始します...\n');
params_A = params_base;
params_A.F_static_max = F_static_required; % ★ 静的荷重を設定
params_A.t_ramp_static = t_ramp_static;   % ★ 静的荷重のランプ時間を設定
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

% --- ケースB: プリテンションあり + 静的荷重あり ---
fprintf('\nケースB（プリテンションあり）のシミュレーションを開始します...\n');
params_B = params_base;
membrane_Area = params_B.w1 * params_B.t1;
strain_y_preload = 2 * z_preload_for_tension^2 / params_B.L^2;
Tension_Force_y = params_B.E1 * strain_y_preload * membrane_Area;
k_tension = 2 * Tension_Force_y / params_B.L;
params_B.k1_effective = params_B.k1 + k_tension;
params_B.F_static_max = F_static_required; % ★ ケースAと同じ静的荷重を設定
params_B.t_ramp_static = t_ramp_static;   % ★ 静的荷重のランプ時間を設定
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

% (fprintf部分は省略)

% -------------------------
% 3. 比較プロット
% -------------------------
figure('Name', 'Frequency Response Comparison');
% (グラフ描画部分は変更なし)
hold on;
plot(freq_range, amp_A*1000, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('ケースA: プリテンションなし (共振点 %.2f Hz)', f_res_A));
plot(freq_range, amp_B*1000, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('ケースB: プリテンションあり (共振点 %.2f Hz)', f_res_B));
plot(f_res_A, max(amp_A)*1000, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
plot(f_res_B, max(amp_B)*1000, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
title('周波数応答の比較'); xlabel('周波数 [Hz]'); ylabel('振幅 [mm]');
legend; grid on; set(gca,'Fontsize',16,'FontWeight','bold');

figure('Name', 'Time Response Comparison');
hold on;
plot(t_A, z_A*1000, 'b-', 'DisplayName', 'ケースA: プリテンションなし');
plot(t_B, z_B*1000, 'r-', 'DisplayName', 'ケースB: プリテンションあり');
xline(t_start_dynamic, 'm--', 'LineWidth', 1.5, 'Label', '音波 開始');
yline(z_static_target*1000, 'k:', 'LineWidth', 1, 'DisplayName', '目標の初期たわみ');
hold off;
title('時間応答の比較 : z方向');
xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
legend; grid on; set(gca,'Fontsize',16,'FontWeight','bold');

% -------------------------
% 4. 関数定義
% -------------------------
function p = generate_params(p)
    p.I1 = p.w1 * p.t1^3 / 12; p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.m1 = p.roh1 * p.L * p.w1 * p.t1; p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1);
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3);
    p.m2 = p.S2 * p.w2 * p.roh2; p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2);
end

% ★ freq_sweepを修正：振幅評価を自動中心検出に戻す
function amp = freq_sweep(p, freq_range, tspan, x0, ode_func_handle)
    amp = zeros(size(freq_range));
    h_waitbar = waitbar(0, '周波数スイープを実行中...');
    for i = 1:length(freq_range)
        waitbar(i/length(freq_range), h_waitbar, sprintf('周波数: %d Hz', freq_range(i)));
        omega = 2*pi*freq_range(i);
        [t_sweep,x_sweep] = ode45(@(t,x) ode_func_handle(t,x,p,omega), tspan, x0);
        z_sweep = x_sweep(:,1);
        final_segment_indices = find(t_sweep > tspan(2)*0.75);
        if isempty(final_segment_indices)
            final_segment = z_sweep(end-100:end);
        else
            final_segment = z_sweep(final_segment_indices);
        end
        center_z = mean(final_segment); % 終盤のデータから振動中心を自動で検出
        amp(i) = max(abs(final_segment - center_z)); % その中心からの最大振幅を評価
    end
    close(h_waitbar);
end

% ケースA用：元の（プリテンションなしの）運動方程式
function dx = nonlinearForcedODE_no_tension(t,x,p,omega)
    z = x(1); dz = x(2);
    ramp_static_factor = min(1.0, t / p.t_ramp_static);
    F_static_current = p.F_static_max * ramp_static_factor; % ★ 静的荷重を追加
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2* p.t2) * p.w1;
        F_dynamic = amplitude * cos(omega * time_since_start);
    end
    F_total = F_dynamic + F_static_current; % ★
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1*z - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end

% ケースB用：プリテンションあり（剛性が高い）の運動方程式
function dx = nonlinearForcedODE_with_tension(t,x,p,omega)
    z = x(1); dz = x(2);
    ramp_static_factor = min(1.0, t / p.t_ramp_static);
    F_static_current = p.F_static_max * ramp_static_factor; % ★ 静的荷重を追加
    F_dynamic = 0;
    if t >= p.t_start_dynamic
        time_since_start = t - p.t_start_dynamic;
        ramp_dynamic_factor = min(1.0, time_since_start / p.t_ramp_dynamic);
        amplitude = ramp_dynamic_factor * p.P0 * (p.L - 2* p.t2) * p.w1;
        F_dynamic = amplitude * cos(omega * time_since_start);
    end
    F_total = F_dynamic + F_static_current; % ★
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ( - (p.c1 + 8*p.c2*z^2/denom)*dz - p.k1_effective*z - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
    dx = [dz; ddz];
end