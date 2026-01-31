% =========================================================================
% 2自由度 シミュレーション (自立維持 + 剛性一定 + 適応的周波数スイープ)【並列化版】
% ・初期たわみ: フレームの内力(y_natural)で自立維持（真の平衡点）
% ・剛性変化: なし (k1は常に一定)
% ・並列計算ツールボックスを使用した高速化
% ・適応的周波数スイープ: 粗いスキャン→ピーク検出→細かいスキャン
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
% --- スイープ設定 ---
freq_coarse = 1:1:100;      % [Hz] 粗いスキャン (1 Hz刻み)
freq_fine_range = 3;        % [Hz] ピーク周辺のスキャン範囲
freq_fine_step = 0.1;       % [Hz] 細かいスキャンの刻み幅
tspan_sweep  = [0 2.0];    % [s]

% --- 音波設定 ---
P0_drive     = 2;         % [Pa]

% --- 本番シミュレーション設定 ---
tspan_main   = [0 5.0];   % [s]
t_sound_on   = 1;       % [s]

% --- 初期たわみ設定 (実験値準拠) ---
x_push = 0.261e-3;          % [m] フレーム押し込み量
C_eff  = 0.2794;          % 有効形状定数

z_target_initial_mm = sqrt(2 * (x_push * 1000) / C_eff);  % [mm]
z_target_initial    = z_target_initial_mm / 1000;         % [m]

% --- 構造・材料パラメータ ---
params = struct( ...
    'L',   18e-3, ...      % [m]
    'w1',  22.5e-3, ...      % [m]
    'w1_stiff', 5.0e-3, ...     % [m] ★剛性計算用（一番細い部分 5mm）
    'w1_mass',  112.5e-3, ...    % [m] ★質量・面積計算用（平均値 22.5mm）
    't1',  0.1e-3, ...    % [m]
    'E1',  3.45e9, ...      % [Pa]
    'rho1',1250, ...       % [kg/m^3]
    'S2',  152.522e-6, ...  % [m^2]
    'rho2',1250, ...       % [kg/m^3]
    'k2_base',222.0, ...   % [N/m]
    'delta',0.05, ...      % 減衰比
    'P0',  P0_drive, ...   % [Pa]
    'freq',0 ...           % [Hz]（後で上書き）
);

% -------------------------------------------------------------------------
% 2. 物理定数と「真の平衡点」の計算
% -------------------------------------------------------------------------
params.m1 = params.rho1 * params.L * params.w1_mass * params.t1;
params.m2 = params.rho2 * params.S2 * params.L;
params.I1 = params.w1_stiff * params.t1^3 / 12;
params.k1       = 384 * params.E1 * params.I1 / (params.L^3);
params.K_couple = params.E1 * (params.w1_stiff * params.t1) / params.L;
params.c1       = 2 * params.delta * sqrt(params.m1 * params.k1);
params.c2       = 2 * params.delta * sqrt(params.m2 * params.k2_base);
params.Gamma = (C_eff / 2) * 1000;   % [1/m]

z_eq = z_target_initial;
coupling_eq = params.k1 / (2 * params.Gamma * params.K_couple);
y_eq = params.Gamma * z_eq^2 + coupling_eq;
y_natural = y_eq + (params.K_couple * coupling_eq) / params.k2_base;

params.z_eq      = z_eq;
params.y_eq      = y_eq;
params.y_natural = y_natural;

fprintf('--- 初期条件 (自立モデル: 真の平衡点) ---\n');
fprintf('フレーム押し込み量:    %.3f mm\n', x_push * 1000);
fprintf('目標初期たわみ z_eq:    %.3f mm\n', z_eq * 1000);
fprintf('対応する平衡点 y_eq:    %.3f mm\n', y_eq * 1000);
fprintf('フレーム自然長 y_nat:   %.3f mm\n', y_natural * 1000);
fprintf('--------------------------------------\n');

x0 = [z_eq; y_eq; 0; 0];

% ODE設定（スイープ用は精度を少し緩める）
options_sweep = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);

% -------------------------------------------------------------------------
% 3. 適応的周波数スイープ
% -------------------------------------------------------------------------
% 並列プールの起動
if isempty(gcp('nocreate'))
    parpool;
end

% === Phase 1: 粗いスキャン (1 Hz刻み) ===
fprintf('Phase 1: 粗いスキャン (1-%d Hz, 1 Hz刻み)...\n', max(freq_coarse));
amp_coarse = zeros(size(freq_coarse));

D1 = parallel.pool.DataQueue;
h_wait1 = waitbar(0, 'Phase 1: 粗いスキャン中...');
afterEach(D1, @(~) waitbar(sum(amp_coarse > 0)/length(freq_coarse), h_wait1));

parfor i = 1:length(freq_coarse)
    f_curr = freq_coarse(i);
    
    params_sweep      = params;
    params_sweep.freq = f_curr;
    
    [~, x_sw] = ode45(@(t,x) equations_2DOF_SelfSustain(t, x, params_sweep, 0), ...
                      tspan_sweep, x0, options_sweep);
    
    z_sw    = x_sw(:,1);
    z_steady = z_sw(round(end/2):end); 
    amp_coarse(i) = (max(z_steady) - min(z_steady)) / 2;
    
    send(D1, i);
end
close(h_wait1);

% ピーク検出
[max_amp_coarse, idx_peak] = max(amp_coarse);
f_peak_coarse = freq_coarse(idx_peak);
fprintf('Phase 1 完了。ピーク検出: %.1f Hz (振幅: %.3f mm)\n', ...
        f_peak_coarse, max_amp_coarse*1000);

% === Phase 2: ピーク周辺の細かいスキャン (0.1 Hz刻み) ===
f_fine_start = max(1, f_peak_coarse - freq_fine_range);
f_fine_end   = min(100, f_peak_coarse + freq_fine_range);
freq_fine = f_fine_start:freq_fine_step:f_fine_end;

fprintf('Phase 2: 細かいスキャン (%.1f-%.1f Hz, %.2f Hz刻み)...\n', ...
        f_fine_start, f_fine_end, freq_fine_step);
amp_fine = zeros(size(freq_fine));

D2 = parallel.pool.DataQueue;
h_wait2 = waitbar(0, 'Phase 2: 細かいスキャン中...');
afterEach(D2, @(~) waitbar(sum(amp_fine > 0)/length(freq_fine), h_wait2));

parfor i = 1:length(freq_fine)
    f_curr = freq_fine(i);
    
    params_sweep      = params;
    params_sweep.freq = f_curr;
    
    [~, x_sw] = ode45(@(t,x) equations_2DOF_SelfSustain(t, x, params_sweep, 0), ...
                      tspan_sweep, x0, options_sweep);
    
    z_sw    = x_sw(:,1);
    z_steady = z_sw(round(end/2):end); 
    amp_fine(i) = (max(z_steady) - min(z_steady)) / 2;
    
    send(D2, i);
end
close(h_wait2);

% 最終的な共振周波数決定
[max_amp, idx_res] = max(amp_fine);
f_res = freq_fine(idx_res);
fprintf('Phase 2 完了。共振周波数: %.2f Hz (振幅: %.3f mm)\n', ...
        f_res, max_amp*1000);

% -------------------------------------------------------------------------
% 4. 本番シミュレーション (共振周波数にて)
% -------------------------------------------------------------------------
fprintf('共振周波数での時間応答を計算中...\n');
params.freq = f_res; 
options     = odeset('RelTol',1e-6, 'AbsTol',1e-9);
[t, x] = ode45(@(t,x) equations_2DOF_SelfSustain(t, x, params, t_sound_on), ...
               tspan_main, x0, options);

z = x(:,1);
y = x(:,2);

% -------------------------------------------------------------------------
% 5. プロット
% -------------------------------------------------------------------------
% --- Figure 1: 周波数応答（2段階スキャン結果） ---
figure('Position', [50, 100, 800, 400], 'Color', 'w');

% 粗いスキャンと細かいスキャンを両方表示
plot(freq_coarse, amp_coarse*1000, 'b.-', 'LineWidth', 1.0, 'MarkerSize', 6);
hold on;
plot(freq_fine, amp_fine*1000, 'r.-', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(f_res, max_amp*1000, 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'g');
hold off;

title('適応的周波数スイープ (剛性一定・自立モデル)');
xlabel('周波数 [Hz]'); ylabel('Z方向 振幅 [mm]');
grid on; 
legend('Phase 1: 粗いスキャン (1 Hz)', ...
       sprintf('Phase 2: 細かいスキャン (%.2f Hz)', freq_fine_step), ...
       sprintf('共振点: %.2f Hz', f_res), ...
       'Location', 'best');

% --- Figure 2: 時間応答 ---
figure('Position', [700, 100, 800, 600], 'Color', 'w');

% (1) Z方向 (膜)
subplot(2,1,1);
plot(t, z*1000, 'b-');
ylabel('Z 変位 [mm]'); 
title(sprintf('膜のたわみ (Z)(駆動: %.2f Hz)', f_res)); 
grid on;

% (2) Y方向 (フレーム)
subplot(2,1,2);
y_extension = (params.y_eq - y) * 1000;
plot(t, y_extension, 'r-', 'LineWidth', 1.5);
ylabel('拡張量 [mm]'); 
title('フレームの変位 (Y)'); 
grid on;

% =========================================================================
% 2自由度 運動方程式 (自立維持・剛性一定・音波あり)
% =========================================================================
function dx = equations_2DOF_SelfSustain(t, x, p, t_sound_on)
    z  = x(1); 
    y  = x(2); 
    dz = x(3); 
    dy = x(4);
    
    k1_curr = p.k1;
    
    coupling_term = y - (p.Gamma * z^2);
    
    F_restore_z = k1_curr * z ...
                  - p.K_couple * coupling_term * (2 * p.Gamma * z);
    F_restore_y = p.k2_base * (y - p.y_natural) ...
                  + p.K_couple * coupling_term;
    
    if t >= t_sound_on
        ramp    = min(1.0, (t - t_sound_on)/0.1);
        F_sound = ramp * p.P0 * (p.L * p.w1_mass) * cos(2 * pi * p.freq * t);
    else
        F_sound = 0;
    end
    
    ddz = (F_sound - p.c1 * dz - F_restore_z) / p.m1;
    ddy = (0 - p.c2 * dy - F_restore_y)     / p.m2;
    
    dx = [dz; dy; ddz; ddy];
end
