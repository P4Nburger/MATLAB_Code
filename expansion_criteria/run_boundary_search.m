%RUN_BOUNDARY_SEARCH  スナップスルーが起きる押し込み量の境界を自動で探索する
%
%   押し込み量 x_push を変えながら
%     「周波数スイープで共振周波数を求める → その周波数で時間応答を計算する」
%   を繰り返し、スナップスルーが起きる押し込み量の境界を求める。
%
%     Phase 0  : 二分探索で境界を 0.05 mm の幅まで絞る
%                （高速化のため、周波数スイープの2段目は省略する）
%     Phase 1〜: 刻み幅を 0.01 → 0.001 → 0.0001 mm と細かくして境界を詰める
%                （境界が 0.1 mm 以上なら 0.001 mm、未満なら 0.0001 mm で終了）
%
%   スナップスルーの判定: 時間応答で min z(t) < -0.8 * z_eq
%
%   旧ファイル名: run_simulation_optimized.m
%   論文: 第4章 図4.2
%
%   See also RUN_RESPONSE_FAST.

clear; clc; close all;

%% 1. 設定
p = model_params('P0', 2);    % 音圧振幅 2 Pa。その他の値は model_params.m を参照

% 周波数スイープ
sweep = struct( ...
    'freqs',        1:1:100, ...  % [Hz] 1段目（1 Hz 刻み）
    'refine_range', 3, ...        % [Hz] 2段目はピークの前後 ±3 Hz
    'refine_step',  0.1, ...      % [Hz] 2段目の刻み幅
    'tspan',        [0 2.0], ...  % [s]  1周波数あたりの計算時間
    'ode_opts',     odeset('RelTol', 1e-4, 'AbsTol', 1e-7), ...
    'use_parallel', start_parallel_pool());

% 時間応答とスナップスルーの判定
sim = struct( ...
    'tspan',      [0 5.0], ...    % [s] 計算時間
    't_sound_on', 1, ...          % [s] 音波を入れる時刻
    'ode_opts',   odeset('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'threshold',  0.8);           % min z < -threshold * z_eq ならスナップスルーあり

% 探索範囲（必要に応じて広げる）
x_search_min    = 0.02e-3;   % [m] スナップスルーが確実に起きる押し込み量
x_search_max    = 0.30e-3;   % [m] スナップスルーが確実に起きない押し込み量
rough_precision = 0.05e-3;   % [m] Phase 0 で絞り込む幅

%% 2. Phase 0: 二分探索で探索範囲を決める
fprintf('========== Phase 0: 探索範囲の決定 (高速版) ==========\n');

fprintf('上限チェック: %.3f mm ... ', x_search_max * 1000);
r = check_snapthrough(x_search_max, p, sweep, sim, true);
if r.is_snap
    error('設定した上限 (%.3f mm) でもスナップしました。探索範囲の上限を上げてください。', x_search_max * 1000);
end
fprintf('なし (OK)\n');

fprintf('下限チェック: %.3f mm ... ', x_search_min * 1000);
r = check_snapthrough(x_search_min, p, sweep, sim, true);
if ~r.is_snap
    error('設定した下限 (%.3f mm) でもスナップしませんでした。探索範囲の下限を下げてください。', x_search_min * 1000);
end
fprintf('あり (OK)\n');

low  = x_search_min;
high = x_search_max;
while (high - low) > rough_precision
    mid = (low + high) / 2;
    fprintf('二分探索中: %.3f mm ... ', mid * 1000);
    r = check_snapthrough(mid, p, sweep, sim, true);
    if r.is_snap
        fprintf('あり -> 下限引き上げ\n');
        low = mid;
    else
        fprintf('なし -> 上限引き下げ\n');
        high = mid;
    end
end

fprintf('探索初期範囲決定: %.3f mm (あり) ～ %.3f mm (なし)\n', low * 1000, high * 1000);

%% 3. Phase 1〜: 刻み幅を細かくしながら境界を詰める
fprintf('\n========== 多段階詳細探索開始 ==========\n');

current_lower   = low;     % [m]  スナップスルーあり
current_upper   = high;    % [m]  スナップスルーなし
current_step_mm = 0.05;    % [mm] 現在の精度
phase_count     = 1;

while true
    % 境界の大きさに応じて目標精度を決める
    if current_lower * 1000 >= 0.1 - 1e-9
        required_step = 0.001;
        order_msg = '>= 0.1 mm (目標 0.001)';
    else
        required_step = 0.0001;
        order_msg = '< 0.1 mm (目標 0.0001)';
    end

    if current_step_mm <= required_step + 1e-9
        fprintf('>> 現在の精度(%.4f)は、値の大きさ(%s)に対して十分です。探索を終了します。\n', ...
                current_step_mm, order_msg);
        break;
    end

    % 次の刻み幅（0.05 → 0.01 → 0.001 → 0.0001 mm）
    next_step_mm = current_step_mm / 10;
    if next_step_mm < required_step
        next_step_mm = required_step;
    end
    if current_step_mm > 0.01 + 1e-9
        next_step_mm = 0.01;
    end

    fprintf('\n--- Phase %d: %.4f mm 刻みで探索 ---\n', phase_count, next_step_mm);

    start_int = round(current_lower * 1000 / next_step_mm);
    end_int   = round(current_upper * 1000 / next_step_mm);
    x_grid_mm = (start_int:1:end_int) * next_step_mm;
    x_grid    = x_grid_mm * 1e-3;

    fprintf('範囲: %.4f mm ～ %.4f mm (最大%d点)\n', x_grid_mm(1), x_grid_mm(end), length(x_grid));

    % 下から順に調べ、初めて「なし」になったところを新しい境界にする
    found_boundary = false;
    new_lower = [];
    new_upper = [];
    history_x    = [];
    history_snap = [];

    for i = 1:length(x_grid)
        x_curr = x_grid(i);

        if abs(x_curr - current_lower) < 1e-9
            is_snap = true;     % 両端は結果がわかっているので計算しない
        elseif abs(x_curr - current_upper) < 1e-9
            is_snap = false;
        else
            fprintf('  Check: %.5f mm ... ', x_curr * 1000);
            r = check_snapthrough(x_curr, p, sweep, sim, false);
            is_snap = r.is_snap;
            if is_snap, fprintf('あり\n'); else, fprintf('なし -> 境界！\n'); end
        end

        history_x    = [history_x; x_curr];      %#ok<AGROW>
        history_snap = [history_snap; is_snap];  %#ok<AGROW>

        if is_snap
            new_lower = x_curr;
        else
            new_upper = x_curr;
            found_boundary = true;
            break;
        end
    end

    if ~found_boundary
        error('境界が見つかりませんでした。');
    end

    current_lower   = new_lower;
    current_upper   = new_upper;
    current_step_mm = next_step_mm;
    phase_count     = phase_count + 1;
end

%% 4. 最終結果
x_final_snap   = current_lower;
x_final_nosnap = current_upper;
step_final     = current_step_mm;

fprintf('\n========== 最終結果 ==========\n');
fprintf('境界(あり): %.5f mm\n', x_final_snap * 1000);
fprintf('境界(なし): %.5f mm\n', x_final_nosnap * 1000);
fprintf('精度:      %.5f mm\n', step_final);

%% 5. プロット
fprintf('\nプロット生成中...\n');
r_snap   = check_snapthrough(x_final_snap,   p, sweep, sim, false);
r_nosnap = check_snapthrough(x_final_nosnap, p, sweep, sim, false);

figure('Position', [100, 100, 1000, 700], 'Color', 'w');

% (1) 最後の Phase の探索履歴（赤: あり、青: なし）
subplot(2, 2, 1);
colors_plot = repmat([0 0 1], length(history_snap), 1);
colors_plot(logical(history_snap), :) = repmat([1 0 0], nnz(history_snap), 1);
scatter(history_x * 1000, ones(size(history_x)), 100, colors_plot, 'filled');
hold on;
plot([x_final_snap,   x_final_snap]   * 1000, [0.95, 1.05], 'r-', 'LineWidth', 2);
plot([x_final_nosnap, x_final_nosnap] * 1000, [0.95, 1.05], 'b-', 'LineWidth', 2);
hold off;
ylim([0.9 1.1]);
xlabel('初期押込量 [mm]');
title(sprintf('探索履歴 (精度: %.4f mm)', step_final));
grid on;

% (2) 境界（あり）の時間応答
subplot(2, 2, 2);
plot_time_response(r_snap, sim.threshold, 'r');
title(sprintf('あり: %.5f mm @ %.2f Hz', x_final_snap * 1000, r_snap.f_res));

% (3) 境界（なし）の時間応答
subplot(2, 2, 3);
plot_time_response(r_nosnap, sim.threshold, 'b');
title(sprintf('なし: %.5f mm @ %.2f Hz', x_final_nosnap * 1000, r_nosnap.f_res));

% (4) 両者の比較
subplot(2, 2, 4);
plot(r_snap.t, r_snap.z * 1000, 'r-', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('あり (%.5f)', x_final_snap * 1000));
hold on;
plot(r_nosnap.t, r_nosnap.z * 1000, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('なし (%.5f)', x_final_nosnap * 1000));
yline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off;
xlabel('時間 [s]'); ylabel('Z変位 [mm]');
title('境界値での比較');
legend('Location', 'best');
grid on;

%% ローカル関数
function r = check_snapthrough(x_push, p, sweep, sim, fast_mode)
%CHECK_SNAPTHROUGH  押し込み量 x_push でスナップスルーが起きるかを調べる
%   fast_mode = true のときは周波数スイープの2段目を省略する。
[p, x0] = initial_state(p, x_push);
if fast_mode
    sweep.refine_range = [];
end
res = find_resonance(p, x0, sweep);
[t, x] = simulate_response(p, x0, res.f_res, sim.tspan, sim.t_sound_on, sim.ode_opts);

r.t       = t;
r.z       = x(:, 1);
r.z_eq    = p.z_eq;
r.f_res   = res.f_res;
r.is_snap = (min(r.z) < -p.z_eq * sim.threshold);
end

function plot_time_response(r, threshold, color)
%PLOT_TIME_RESPONSE  膜のたわみの時間応答を、反転位置と判定閾値の線と一緒に描く
plot(r.t, r.z * 1000, [color '-'], 'LineWidth', 1.5);
hold on;
yline(0, 'k--', 'LineWidth', 0.5);
yline(-r.z_eq * 1000, 'g:', 'LineWidth', 1.5, 'Label', '-z_{eq}');
yline(-threshold * r.z_eq * 1000, 'm--', 'LineWidth', 1.0, 'Label', '判定閾値');
hold off;
xlabel('時間 [s]'); ylabel('Z変位 [mm]');
grid on;
end
