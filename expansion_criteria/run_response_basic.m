%RUN_RESPONSE_BASIC  1条件での周波数スイープと時間応答（基本版）
%
%   押し込み量を1つ決めて、
%     1) 周波数スイープで共振周波数を求め、
%     2) その周波数の音波を当てたときの時間応答を計算する。
%   1段階のスイープを逐次計算で行う基本版。
%   広い周波数範囲を細かく調べるときは run_response_fast.m を使う。
%
%   旧ファイル名: gouseihenkou.m
%   論文: 第4章
%
%   See also RUN_RESPONSE_FAST, RUN_BOUNDARY_SEARCH.

clear; clc; close all;

%% 1. 設定
% 押し込み量（12/4 時点: 0.07 mm でスナップスルーあり、0.08 mm でなし）
x_push = 0.07e-3;             % [m] フレームの押し込み量

p = model_params('P0', 2);    % 音圧振幅 2 Pa。その他の値は model_params.m を参照

% 周波数スイープ
sweep = struct( ...
    'freqs',        20:1:40, ...  % [Hz]
    'refine_range', [], ...       % 2段目は行わない
    'refine_step',  [], ...
    'tspan',        [0 2.0], ...  % [s] 1周波数あたりの計算時間
    'ode_opts',     [], ...       % ode45 の既定の精度
    'use_parallel', false);

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

%% 3. 周波数スイープ
fprintf('周波数スイープを実行中...\n');
res   = find_resonance(p, x0, sweep);
f_res = res.f_res;
fprintf('スイープ完了。共振周波数: %.2f Hz\n', f_res);

%% 4. 共振周波数での時間応答
fprintf('共振周波数での時間応答を計算中...\n');
[t, x] = simulate_response(p, x0, f_res, tspan_main, t_sound_on, ode_opts_main);
z = x(:, 1);
y = x(:, 2);

%% 5. プロット
% 周波数応答
figure('Position', [50, 100, 600, 400], 'Color', 'w');
plot(res.freqs1, res.amp1 * 1000, 'b.-', 'LineWidth', 1.5);
hold on;
plot(f_res, res.amp_res * 1000, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
title('周波数応答 (剛性一定・自立モデル)');
xlabel('周波数 [Hz]'); ylabel('Z方向 振幅 [mm]');
grid on; legend('応答曲線', sprintf('共振点: %.1f Hz', f_res));

% 時間応答
figure('Position', [700, 100, 800, 600], 'Color', 'w');

subplot(2, 1, 1);   % 膜のたわみ
plot(t, z * 1000, 'b-');
xlabel('時間 [s]'); ylabel('Z 変位 [mm]');
title(sprintf('膜のたわみ (Z)(駆動: %.2f Hz)', f_res));
grid on; xlim([0 5]); ylim([-1.2 1.2]);

subplot(2, 1, 2);   % フレームの変位（拡張量 = 釣り合い位置からの縮み）
y_extension = (p.y_eq - y) * 1000;
plot(t, y_extension, 'r-', 'LineWidth', 1.5);
xlabel('時間 [s]'); ylabel('拡張量 [mm]');
title('フレームの変位 (Y)');
grid on; xlim([0 5]); ylim([-0.1 0.1]);
