%RUN_RESPONSE_FAST  1条件での周波数スイープと時間応答（高速版）
%
%   run_response_basic.m と同じ解析を、
%     ・2段階の周波数スイープ（1 Hz 刻みで広く → ピークの前後を細かく）
%     ・並列計算（Parallel Computing Toolbox があれば）
%   で行う。
%
%   旧ファイル名: scan.m
%   論文: 第4章
%
%   卒論提出後の音圧 5 Pa の検討（旧 onatuhenkou.m）は、
%   P0 = 5, x_push = 0.0735e-3, refine_step = 0.1 にすると再現できる。
%
%   See also RUN_RESPONSE_BASIC, RUN_BOUNDARY_SEARCH.

clear; clc; close all;

%% 1. 設定
x_push = 0.073e-3;            % [m] フレームの押し込み量

p = model_params('P0', 2);    % 音圧振幅 2 Pa。その他の値は model_params.m を参照

% 周波数スイープ
sweep = struct( ...
    'freqs',        1:1:100, ...  % [Hz] 1段目（1 Hz 刻み）
    'refine_range', 3, ...        % [Hz] 2段目はピークの前後 ±3 Hz
    'refine_step',  0.01, ...     % [Hz] 2段目の刻み幅
    'tspan',        [0 2.0], ...  % [s]  1周波数あたりの計算時間
    'ode_opts',     odeset('RelTol', 1e-4, 'AbsTol', 1e-7), ...  % スイープ用に精度を少し緩める
    'use_parallel', start_parallel_pool());

% 時間応答
tspan_main    = [0 5.0];          % [s] 計算時間
t_sound_on    = 1;                % [s] 音波を入れる時刻
ode_opts_main = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

%% 2. 初期状態（釣り合い点）
[p, x0] = initial_state(p, x_push);

fprintf('--- 初期条件 (自立モデル: 真の平衡点) ---\n');
fprintf('フレーム押し込み量:    %.3f mm\n', x_push * 1000);
fprintf('目標初期たわみ z_eq:    %.3f mm\n', p.z_eq * 1000);
fprintf('対応する平衡点 y_eq:    %.3f mm\n', p.y_eq * 1000);
fprintf('フレーム自然長 y_nat:   %.3f mm\n', p.y_natural * 1000);
fprintf('--------------------------------------\n');

%% 3. 2段階の周波数スイープ
fprintf('周波数スイープを実行中...\n');
res   = find_resonance(p, x0, sweep);
f_res = res.f_res;
fprintf('1段目: %g〜%g Hz → ピーク %.1f Hz (振幅: %.3f mm)\n', ...
        res.freqs1(1), res.freqs1(end), res.f_peak1, res.amp_peak1 * 1000);
fprintf('2段目: %.2f〜%.2f Hz (%.2f Hz 刻み) → 共振周波数 %.2f Hz (振幅: %.3f mm)\n', ...
        res.freqs2(1), res.freqs2(end), sweep.refine_step, f_res, res.amp_res * 1000);

%% 4. 共振周波数での時間応答
fprintf('共振周波数での時間応答を計算中...\n');
[t, x] = simulate_response(p, x0, f_res, tspan_main, t_sound_on, ode_opts_main);
z = x(:, 1);
y = x(:, 2);

%% 5. プロット
% 周波数応答（2段階スイープの結果）
figure('Position', [50, 100, 800, 400], 'Color', 'w');
plot(res.freqs1, res.amp1 * 1000, 'b.-', 'LineWidth', 1.0, 'MarkerSize', 6);
hold on;
plot(res.freqs2, res.amp2 * 1000, 'r.-', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(f_res, res.amp_res * 1000, 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'g');
hold off;
title('適応的周波数スイープ (剛性一定・自立モデル)');
xlabel('周波数 [Hz]'); ylabel('Z方向 振幅 [mm]');
grid on;
legend('1段目: 粗いスキャン (1 Hz)', ...
       sprintf('2段目: 細かいスキャン (%.2f Hz)', sweep.refine_step), ...
       sprintf('共振点: %.2f Hz', f_res), ...
       'Location', 'best');

% 時間応答
figure('Position', [700, 100, 800, 600], 'Color', 'w');

subplot(2, 1, 1);   % 膜のたわみ
plot(t, z * 1000, 'b-');
xlabel('時間 [s]'); ylabel('Z 変位 [mm]');
title(sprintf('膜のたわみ (Z)(駆動: %.2f Hz)', f_res));
grid on;

subplot(2, 1, 2);   % フレームの変位（拡張量 = 釣り合い位置からの縮み）
y_extension = (p.y_eq - y) * 1000;
plot(t, y_extension, 'r-', 'LineWidth', 1.5);
xlabel('時間 [s]'); ylabel('拡張量 [mm]');
title('フレームの変位 (Y)');
grid on;
