% =========================================================================
% スナップスルー境界値自動探索プログラム（プロット復刻版）
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. パラメータ設定
% -------------------------------------------------------------------------
C_eff  = 0.2794;
snapthrough_threshold = 0.8;

% --- シミュレーション設定 ---
freq_coarse = 1:1:100;
freq_fine_range = 3;
freq_fine_step = 0.1;
tspan_sweep  = [0 2.0];
P0_drive = 2;
tspan_main = [0 5.0];
t_sound_on = 1;

% --- 構造パラメータ（ユーザー入力値を維持） ---
params = struct('L',18e-3, 'w1',22.5e-3, 'w1_stiff',5.0e-3, 'w1_mass',11.25e-3, ...
    't1',0.1e-3, 'E1',3.45e9, 'rho1',1250, 'S2',152.522e-6, 'rho2',1250, ...
    'k2_base',222.0, 'delta',0.05, 'P0',P0_drive, 'freq',0);

% 物理定数
params.m1 = params.rho1 * params.L * params.w1_mass * params.t1;
params.m2 = params.rho2 * params.S2 * params.L;
params.I1 = params.w1_stiff * params.t1^3 / 12;
params.k1 = 384 * params.E1 * params.I1 / (params.L^3);
params.K_couple = params.E1 * (params.w1_stiff * params.t1) / params.L;
params.c1 = 2 * params.delta * sqrt(params.m1 * params.k1);
params.c2 = 2 * params.delta * sqrt(params.m2 * params.k2_base);
params.Gamma = (C_eff / 2) * 1000;

options_sweep = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);
options_main = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

if isempty(gcp('nocreate')), parpool; end

% -------------------------------------------------------------------------
% 2. Phase 0: 探索範囲の自動決定
% -------------------------------------------------------------------------
fprintf('========== Phase 0: 探索範囲の決定 ==========\n');
x_test_values = (0.02:0.05:0.25) * 1e-3; 
x_lower_limit = []; x_upper_limit = [];

for x_val = x_test_values
    fprintf('テスト: %.3f mm ... ', x_val*1000);
    is_snap = check_snapthrough(x_val, params, C_eff, snapthrough_threshold, freq_coarse, freq_fine_range, freq_fine_step, tspan_sweep, tspan_main, t_sound_on, options_sweep, options_main);
    if is_snap
        x_lower_limit = x_val; fprintf('あり\n');
    else
        x_upper_limit = x_val; fprintf('なし (範囲確保)\n'); break;
    end
end

if isempty(x_lower_limit) || isempty(x_upper_limit)
    if isempty(x_lower_limit)
        fprintf('下限探索中...\n');
        x_back = x_upper_limit - 0.05e-3;
        while x_back > 0
             fprintf('テスト: %.3f mm ... ', x_back*1000);
             if check_snapthrough(x_back, params, C_eff, snapthrough_threshold, freq_coarse, freq_fine_range, freq_fine_step, tspan_sweep, tspan_main, t_sound_on, options_sweep, options_main)
                 x_lower_limit = x_back; fprintf('あり (下限確保)\n'); break;
             else
                 fprintf('なし\n'); x_back = x_back - 0.05e-3;
             end
        end
    end
end
if isempty(x_lower_limit), error('有効な探索範囲が見つかりませんでした。'); end

fprintf('探索初期範囲: %.3f mm (あり) ～ %.3f mm (なし)\n', x_lower_limit*1000, x_upper_limit*1000);

% -------------------------------------------------------------------------
% 3. 多段階詳細探索（動的終了判定）
% -------------------------------------------------------------------------
fprintf('\n========== 多段階詳細探索開始 ==========\n');

current_lower = x_lower_limit;
current_upper = x_upper_limit;
current_step_mm = 0.05; 
phase_count = 1;

while true
    val_mm = current_lower * 1000;
    
    % 動的精度判定
    if val_mm >= 0.1 - 1e-9
        required_step = 0.001; 
        order_msg = '>= 0.1 mm (目標 0.001)';
    else
        required_step = 0.0001; 
        order_msg = '< 0.1 mm (目標 0.0001)';
    end
    
    if current_step_mm <= required_step + 1e-9
        fprintf('>> 現在の精度(%.4f)は、値の大きさ(%s)に対して十分です。探索を終了します。\n', current_step_mm, order_msg);
        break;
    end
    
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
    x_grid_mm = (start_int : 1 : end_int) * next_step_mm;
    x_grid = x_grid_mm * 1e-3;
    
    fprintf('範囲: %.4f mm ～ %.4f mm (最大%d点)\n', x_grid_mm(1), x_grid_mm(end), length(x_grid));
    
    found_boundary = false;
    new_lower = [];
    new_upper = [];
    
    history_x = [];
    history_snap = [];
    
    for i = 1:length(x_grid)
        x_curr = x_grid(i);
        
        if abs(x_curr - current_lower) < 1e-9
            is_snap = true; 
        elseif abs(x_curr - current_upper) < 1e-9
            is_snap = false;
        else
            fprintf('  Check: %.5f mm ... ', x_curr*1000);
            is_snap = check_snapthrough(x_curr, params, C_eff, snapthrough_threshold, freq_coarse, freq_fine_range, freq_fine_step, tspan_sweep, tspan_main, t_sound_on, options_sweep, options_main);
            if is_snap, fprintf('あり\n'); else, fprintf('なし -> 境界！\n'); end
        end
        
        history_x = [history_x; x_curr];
        history_snap = [history_snap; is_snap];
        
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
    
    current_lower = new_lower;
    current_upper = new_upper;
    current_step_mm = next_step_mm;
    
    phase_count = phase_count + 1;
end

% -------------------------------------------------------------------------
% 4. 最終結果とプロット (ここを旧バージョンの詳細プロットに戻しました)
% -------------------------------------------------------------------------
x_final_snap = current_lower;
x_final_nosnap = current_upper;
step_final = current_step_mm;

fprintf('\n========== 最終結果 ==========\n');
fprintf('境界(あり): %.5f mm\n', x_final_snap*1000);
fprintf('境界(なし): %.5f mm\n', x_final_nosnap*1000);
fprintf('精度:      %.5f mm\n', step_final);

% プロット用データ生成
fprintf('\nプロット生成中...\n');
[t_snap, z_snap, f_res_snap] = run_sim(x_final_snap, params, C_eff, freq_coarse, freq_fine_range, freq_fine_step, tspan_sweep, tspan_main, t_sound_on, options_sweep, options_main);
[t_nosnap, z_nosnap, f_res_nosnap] = run_sim(x_final_nosnap, params, C_eff, freq_coarse, freq_fine_range, freq_fine_step, tspan_sweep, tspan_main, t_sound_on, options_sweep, options_main);

figure('Position', [100, 100, 1000, 700], 'Color', 'w');

% (1) 詳細探索の結果
subplot(2,2,1);
colors_plot = zeros(length(history_snap), 3);
for i = 1:length(history_snap)
    if history_snap(i)
        colors_plot(i,:) = [1 0 0];
    else
        colors_plot(i,:) = [0 0 1];
    end
end
scatter(history_x*1000, ones(size(history_x)), 100, colors_plot, 'filled');
hold on;
plot([x_final_snap*1000, x_final_snap*1000], [0.95, 1.05], 'r-', 'LineWidth', 2);
plot([x_final_nosnap*1000, x_final_nosnap*1000], [0.95, 1.05], 'b-', 'LineWidth', 2);
hold off;
ylim([0.9 1.1]);
xlabel('初期押込量 [mm]');
title(sprintf('探索履歴 (精度: %.4f mm)', step_final));
grid on;

% (2) 境界値（あり）の時間応答
subplot(2,2,2);
z_eq_lower = sqrt(2 * (x_final_snap * 1000) / C_eff) / 1000;
plot(t_snap, z_snap*1000, 'r-', 'LineWidth', 1.5);
hold on;
yline(0, 'k--', 'LineWidth', 0.5);
yline(-z_eq_lower*1000, 'g:', 'LineWidth', 1.5, 'Label', '閾値');
hold off;
xlabel('時間 [s]'); ylabel('Z変位 [mm]');
title(sprintf('あり: %.5f mm @ %.2f Hz', x_final_snap*1000, f_res_snap));
grid on;

% (3) 境界値（なし）の時間応答
subplot(2,2,3);
z_eq_upper = sqrt(2 * (x_final_nosnap * 1000) / C_eff) / 1000;
plot(t_nosnap, z_nosnap*1000, 'b-', 'LineWidth', 1.5);
hold on;
yline(0, 'k--', 'LineWidth', 0.5);
yline(-z_eq_upper*1000, 'g:', 'LineWidth', 1.5, 'Label', '閾値');
hold off;
xlabel('時間 [s]'); ylabel('Z変位 [mm]');
title(sprintf('なし: %.5f mm @ %.2f Hz', x_final_nosnap*1000, f_res_nosnap));
grid on;

% (4) 両者の比較
subplot(2,2,4);
plot(t_snap, z_snap*1000, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('あり (%.5f)', x_final_snap*1000));
hold on;
plot(t_nosnap, z_nosnap*1000, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('なし (%.5f)', x_final_nosnap*1000));
yline(0, 'k--', 'LineWidth', 0.5);
hold off;
xlabel('時間 [s]'); ylabel('Z変位 [mm]');
title('境界値での比較');
legend('Location', 'best');
grid on;

% =========================================================================
% 関数定義
% =========================================================================
function is_snap = check_snapthrough(x_push, params, C_eff, threshold, f_c, f_fr, f_fs, t_sw, t_m, t_on, opt_sw, opt_m)
    [~, z_out, ~] = run_sim(x_push, params, C_eff, f_c, f_fr, f_fs, t_sw, t_m, t_on, opt_sw, opt_m);
    z_eq = sqrt(2*(x_push*1000)/C_eff)/1000;
    is_snap = (min(z_out) < -z_eq * threshold);
end

% run_simを修正して周波数(f_res)も返すように変更
function [t, z, f_res] = run_sim(x_push, params, C_eff, f_c, f_fr, f_fs, t_sw, t_m, t_on, opt_sw, opt_m)
    z_eq = sqrt(2*(x_push*1000)/C_eff)/1000;
    coup_eq = params.k1/(2*params.Gamma*params.K_couple);
    y_eq = params.Gamma*z_eq^2 + coup_eq;
    y_nat = y_eq + (params.K_couple*coup_eq)/params.k2_base;
    p = params; p.z_eq=z_eq; p.y_eq=y_eq; p.y_natural=y_nat;
    x0=[z_eq; y_eq; 0; 0];
    
    amp_c = zeros(size(f_c));
    parfor j=1:length(f_c), p_s=p; p_s.freq=f_c(j); [~,x]=ode45(@(t,x) eqn(t,x,p_s,0), t_sw, x0, opt_sw); z_s=x(round(end/2):end,1); amp_c(j)=(max(z_s)-min(z_s))/2; end
    [~,idx]=max(amp_c); f_peak=f_c(idx);
    
    f_fine = max(1, f_peak-f_fr):f_fs:min(100, f_peak+f_fr);
    amp_f = zeros(size(f_fine));
    parfor j=1:length(f_fine), p_s=p; p_s.freq=f_fine(j); [~,x]=ode45(@(t,x) eqn(t,x,p_s,0), t_sw, x0, opt_sw); z_s=x(round(end/2):end,1); amp_f(j)=(max(z_s)-min(z_s))/2; end
    [~,idx]=max(amp_f); f_res = f_fine(idx); % 共振周波数を取得
    
    p.freq = f_res;
    [t, x_out] = ode45(@(t,x) eqn(t, x, p, t_on), t_m, x0, opt_m);
    z = x_out(:,1);
end

function dx = eqn(t, x, p, t_on)
    z=x(1); y=x(2); dz=x(3); dy=x(4);
    coup = y - (p.Gamma * z^2);
    F_z = p.k1*z - p.K_couple*coup*(2*p.Gamma*z);
    F_y = p.k2_base*(y-p.y_natural) + p.K_couple*coup;
    F_s=0; if t>=t_on, F_s=min(1,(t-t_on)/0.1)*p.P0*(p.L*p.w1_mass)*cos(2*pi*p.freq*t); end
    dx = [dz; dy; (F_s - p.c1*dz - F_z)/p.m1; (0 - p.c2*dy - F_y)/p.m2];
end