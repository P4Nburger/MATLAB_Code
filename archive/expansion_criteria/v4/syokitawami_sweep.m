% =========================================================================
% 2自由度 拡張シミュレーション (v2改 + 周波数スイープ実装版)
% ・初期たわみ: フレーム押し込み量から算出
% ・周波数スイープ: 共振周波数を探索してプロット
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
% --- スイープ設定 ---
freq_range = 20:1:30;  % [Hz] スイープする周波数範囲
tspan_sweep = [0 2.0];  % [s] スイープ時の各周波数の計算時間（短めでOK）

% --- 本番シミュレーション設定 ---
tspan_main = [0 1.0];   % [s] 共振周波数での詳細シミュレーション時間
t_sound_on = 0.2;       % [s] 音波開始

% --- 初期たわみ設定 (実験値準拠) ---
x_push = 0.1e-3;        % [m] フレーム押し込み量
C_eff = 0.2794;         % 有効形状定数
z_target_initial_mm = sqrt(2 * (x_push * 1000) / C_eff);
z_target_initial = z_target_initial_mm / 1000; % [m]

% --- 音圧設定 ---
P0_drive = 25;          % [Pa]

% --- 剛性変化パラメータ ---
stiffness_factor = 3000;
stiffness_threshold = 0.02;
amplitude_sensitivity = 3;
amplitude_decay = 5.0;

params = struct( ...
    'L', 30e-3, 'w1', 10e-3, 't1', 0.02e-3, 'E1', 197e9, 'rho1', 7930, ...
    'S2', 126.87e-6, 'rho2', 1240, 'k2_base', 177.0, ...
    'delta', 0.05, 'P0', P0_drive, 'freq', 0, ... % freqは後で設定
    'stiffness_factor', stiffness_factor, ...
    'threshold', stiffness_threshold, ...
    'amp_sens', amplitude_sensitivity, 'amp_decay', amplitude_decay ...
);

% -------------------------------------------------------------------------
% 2. 物理定数と平衡点の計算
% -------------------------------------------------------------------------
params.m1 = params.rho1 * params.L * params.w1 * params.t1;
params.m2 = params.rho2 * params.S2 * params.L;
params.I1 = params.w1 * params.t1^3 / 12;
params.k1 = 384 * params.E1 * params.I1 / (params.L^3);
params.K_couple = params.E1 * (params.w1 * params.t1) / params.L;
params.c1 = 2 * params.delta * sqrt(params.m1 * params.k1);
params.c2 = 2 * params.delta * sqrt(params.m2 * params.k2_base);

% --- 平衡点の計算 ---
z = z_target_initial;
term = (params.k1 * z * params.L) / (4 * z * params.K_couple);
y_eq = term + (2 * z^2 / params.L);
term_couple = params.K_couple * (y_eq - 2 * z^2 / params.L);
y_natural = y_eq + term_couple / params.k2_base;
params.y_natural = y_natural;
params.z_eq = z;
params.y_eq = y_eq;

fprintf('--- 初期条件 ---\n');
fprintf('フレーム押し込み量: %.3f mm\n', x_push * 1000);
fprintf('計算された初期たわみ: %.3f mm\n', z * 1000);
fprintf('--------------------------------------\n');

% -------------------------------------------------------------------------
% 3. 周波数スイープの実行
% -------------------------------------------------------------------------
fprintf('周波数スイープを実行中 (%.0f Hz ～ %.0f Hz)...\n', min(freq_range), max(freq_range));

% 初期条件
x0 = [z_target_initial; y_eq; 0; 0; 0];

amp_data = zeros(size(freq_range));
h_wait = waitbar(0, '周波数スイープ中...');

for i = 1:length(freq_range)
    f_curr = freq_range(i);
    waitbar(i/length(freq_range), h_wait, sprintf('解析中: %.0f Hz', f_curr));
    
    % パラメータ更新
    params_sweep = params;
    params_sweep.freq = f_curr;
    
    % シミュレーション実行 (音波は最初からONにする)
    % ※スイープ時は過渡応答を早く終わらせるため t_sound_on=0 とする
    [~, x_sw] = ode45(@(t,x) equations_2DOF_Threshold(t, x, params_sweep, 0), tspan_sweep, x0);
    
    % 振幅の計算 (後半のデータを使用)
    z_sw = x_sw(:,1);
    n_data = length(z_sw);
    z_steady = z_sw(round(n_data*0.7):end); % 後半30%
    amp_data(i) = (max(z_steady) - min(z_steady)) / 2;
end
close(h_wait);

% 共振周波数の特定
[max_amp, idx_res] = max(amp_data);
f_res = freq_range(idx_res);

fprintf('スイープ完了。\n');
fprintf('検出された共振周波数: %.2f Hz (最大振幅: %.3f mm)\n', f_res, max_amp*1000);
fprintf('--------------------------------------\n');

% -------------------------------------------------------------------------
% 4. 本番シミュレーション (共振周波数にて)
% -------------------------------------------------------------------------
fprintf('共振周波数 (%.2f Hz) での時間応答を計算中...\n', f_res);

params.freq = f_res; % 周波数を共振点に設定
options = odeset('RelTol',1e-6, 'AbsTol',1e-9);
[t, x] = ode45(@(t,x) equations_2DOF_Threshold(t, x, params, t_sound_on), tspan_main, x0, options);

z = x(:,1);
y = x(:,2);
A = x(:,5);

% -------------------------------------------------------------------------
% 5. プロット
% -------------------------------------------------------------------------

% --- Figure 1: 周波数応答 ---
figure('Position', [50, 100, 600, 400], 'Color', 'w');
plot(freq_range, amp_data*1000, 'b.-', 'LineWidth', 1.5);
hold on;
plot(f_res, max_amp*1000, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
title('周波数応答 (振幅特性)');
xlabel('周波数 [Hz]'); ylabel('Z方向 振幅 [mm]');
grid on; legend('応答曲線', sprintf('共振点: %.1f Hz', f_res));

% --- Figure 2: 時間応答 (Z, Y, 剛性) ---
figure('Position', [700, 100, 800, 900], 'Color', 'w');

% (1) Z方向 (膜)
subplot(3,1,1);
plot(t, z*1000, 'b-');
hold on; 
xline(t_sound_on, 'r--', 'Label', 'Sound ON');
yline(0, 'k-', 'LineWidth', 1); 
hold off;
ylabel('Z 変位 [mm]'); 
title(sprintf('膜のたわみ (駆動: %.2f Hz)', f_res)); 
grid on;

% (2) Y方向 (フレーム)
subplot(3,1,2);
y_extension = (params.y_eq - y) * 1000;
plot(t, y_extension, 'r-', 'LineWidth', 1.5);
hold on;
xline(t_sound_on, 'r--');
yline(0, 'g:', 'LineWidth', 1.5, 'Label', '初期位置');
hold off;
ylabel('拡張量 [mm]'); 
title('フレームの変位 (Y)'); 
grid on;

% (3) 内部変数 (振動の激しさA と 剛性k1)
subplot(3,1,3);
yyaxis left; 
plot(t, A, 'b-', 'LineWidth', 1.5);
hold on;
yline(params.threshold, 'b:', 'LineWidth', 1.5, 'Label', 'しきい値');
hold off;
ylabel('振動の激しさ A [m/s]');
ylim([0, max(A)*1.2]); 

yyaxis right; 
A_effective = max(0, A - params.threshold);
k1_dynamic = params.k1 .* (1 + params.stiffness_factor * A_effective);
plot(t, k1_dynamic, 'r-', 'LineWidth', 1.5);
ylabel('膜の剛性 k1 [N/m]');
ylim([params.k1*0.9, max(k1_dynamic)*1.1]);
title('振動の激しさと剛性の変化'); 
grid on; xlabel('時間 [s]');


% =========================================================================
% 2自由度 運動方程式 (しきい値付き)
% =========================================================================
function dx = equations_2DOF_Threshold(t, x, p, t_sound_on)
    z  = x(1);
    y  = x(2);
    dz = x(3);
    dy = x(4);
    A  = x(5); 
    
    % 1. 剛性の動的変化
    A_effective = max(0, A - p.threshold);
    k1_curr = p.k1 * (1 + p.stiffness_factor * A_effective);

    % 2. 復元力
    coupling_term = y - (2 * z^2 / p.L);
    
    F_restore_z = k1_curr * z - p.K_couple * coupling_term * (4 * z / p.L);
    F_restore_y = p.k2_base * (y - p.y_natural) + p.K_couple * coupling_term;
    
    % 3. 外力
    if t >= t_sound_on
        F_sound = p.P0 * (p.L * p.w1) * cos(2 * pi * p.freq * t);
    else
        F_sound = 0;
    end
    
    % 4. 運動方程式
    ddz = (F_sound - p.c1 * dz - F_restore_z) / p.m1;
    ddy = (0 - p.c2 * dy - F_restore_y) / p.m2;
    
    % 5. 振幅Aの更新
    dA = p.amp_sens * abs(dz) - p.amp_decay * A;
    
    dx = [dz; dy; ddz; ddy; dA];
end