%% parameter_search_zeq_parallel.m
% 目的：初期たわみ z_eq を最大化（CPU並列計算版）

clear; clc; close all;

%% ========== 並列プール起動 ==========
% 利用可能なコア数を確認
num_cores = feature('numcores');
fprintf('利用可能なCPUコア数: %d\n', num_cores);

% 並列プールを起動（既に起動済みなら再利用）
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', num_cores); % 全コア使用
    fprintf('並列プールを起動しました（ワーカー数: %d）\n', pool.NumWorkers);
else
    fprintf('既存の並列プールを使用します（ワーカー数: %d）\n', pool.NumWorkers);
end

%% ========== 探索範囲設定 ==========
search_space = struct(...
    'x_push_m',    linspace(0.02e-3, 0.20e-3, 8), ...
    'w1_stiff_m',  [0.1 0.5 1 2 5]*1e-3, ...
    'w1_mass_m',   [10 20 40 60 80 100]*1e-3, ...
    'L_m',         [15 18 21 25]*1e-3 ...
);

fixed = struct(...
    't1_m', 0.1e-3, ...
    'E1_Pa', 3.45e9, ...
    'rho1', 1250, ...
    'P0_Pa', 2, ...
    'delta', 0.05, ...
    'C_eff', 0.2794, ...
    'S2_m2', 152.522e-6, ...
    'rho2', 1250, ...
    'k2_base_Npm', 222.0 ...
);

comp_config = struct(...
    'freq_coarse', 10:10:200, ...
    'freq_fine_half', 10, ...
    'tspan_sweep', [0 2.0], ...
    'ode_opt', odeset('RelTol',1e-6,'AbsTol',1e-9) ...
);

%% ========== 探索実行（並列化） ==========
n_designs = prod(structfun(@length, search_space));
fprintf('\n=== 初期たわみ最大化探索（並列計算版）===\n');
fprintf('探索点数: %d\n', n_designs);
fprintf('推定時間: %.1f 分（%dコア並列）\n\n', 2*60/pool.NumWorkers, pool.NumWorkers);

tic;
results = run_grid_search_parallel(search_space, fixed, comp_config);
elapsed = toc;

%% ========== 結果保存・可視化 ==========
save('zeq_search_results_parallel.mat', 'results', 'search_space', 'fixed', 'comp_config');
fprintf('\n結果を zeq_search_results_parallel.mat に保存しました\n');
fprintf('実際の計算時間: %.1f 分（並列効率: %.1f%%）\n', ...
    elapsed/60, 100*(2*60)/(elapsed*pool.NumWorkers));

plot_zeq_only_results(results);

%% ========== 並列化版探索関数 ==========

function results = run_grid_search_parallel(ss, fixed, cfg)
    % 全設計点の組み合わせを事前生成
    [X_push, W_stiff, W_mass, L_span] = ndgrid(...
        ss.x_push_m, ss.w1_stiff_m, ss.w1_mass_m, ss.L_m);
    
    % 1次元配列化（parfor用）
    X_push_vec = X_push(:);
    W_stiff_vec = W_stiff(:);
    W_mass_vec = W_mass(:);
    L_span_vec = L_span(:);
    n_total = length(X_push_vec);
    
    % 結果格納用（preallocate）
    z_eq_arr = zeros(n_total, 1);
    f_res_arr = zeros(n_total, 1);
    A_z_arr = zeros(n_total, 1);
    converged_arr = false(n_total, 1);
    
    % 進捗表示用
    fprintf('並列計算開始...\n');
    progress_interval = ceil(n_total / 20); % 5%ごとに表示
    
    % ========== 並列ループ（parfor） ==========
    parfor i = 1:n_total
        design = struct(...
            'x_push_m', X_push_vec(i), ...
            'w1_stiff_m', W_stiff_vec(i), ...
            'w1_mass_m', W_mass_vec(i), ...
            'L_m', L_span_vec(i), ...
            't1_m', fixed.t1_m, ...
            'E1_Pa', fixed.E1_Pa, ...
            'rho1', fixed.rho1, ...
            'P0_Pa', fixed.P0_Pa, ...
            'delta', fixed.delta, ...
            'C_eff', fixed.C_eff, ...
            'S2_m2', fixed.S2_m2, ...
            'rho2', fixed.rho2, ...
            'k2_base_Npm', fixed.k2_base_Npm ...
        );
        
        try
            out = evaluate_one_design_zeq(design, cfg);
            z_eq_arr(i) = out.z_eq_mm;
            f_res_arr(i) = out.f_res_Hz;
            A_z_arr(i) = out.A_z_mm;
            converged_arr(i) = true;
        catch ME
            warning('設計点 %d で失敗: %s', i, ME.message);
            z_eq_arr(i) = NaN;
            f_res_arr(i) = NaN;
            A_z_arr(i) = NaN;
            converged_arr(i) = false;
        end
        
        % 進捗表示（parfor内ではwaitbar不可なのでfprintf）
        if mod(i, progress_interval) == 0
            fprintf('  進捗: %d/%d (%.0f%%)\n', i, n_total, 100*i/n_total);
        end
    end
    
    % 構造体配列に変換
    results = struct('x_push_m', {}, 'w1_stiff_m', {}, 'w1_mass_m', {}, ...
                     'L_m', {}, 'z_eq_mm', {}, 'f_res_Hz', {}, ...
                     'A_z_mm', {}, 'converged', {});
    
    for i = 1:n_total
        results(i).x_push_m = X_push_vec(i);
        results(i).w1_stiff_m = W_stiff_vec(i);
        results(i).w1_mass_m = W_mass_vec(i);
        results(i).L_m = L_span_vec(i);
        results(i).z_eq_mm = z_eq_arr(i);
        results(i).f_res_Hz = f_res_arr(i);
        results(i).A_z_mm = A_z_arr(i);
        results(i).converged = converged_arr(i);
    end
end

%% ========== 以下、評価関数・運動方程式（変更なし） ==========

function out = evaluate_one_design_zeq(des, cfg)
    % 1) 初期たわみ計算 [file:1]
    z_eq_mm = sqrt(2 * (des.x_push_m*1000) / des.C_eff);
    z_eq = z_eq_mm / 1000;
    
    % 2) 等価パラメータ計算 [file:1]
    p = struct();
    p.L = des.L_m;
    p.w1_stiff = des.w1_stiff_m;
    p.w1_mass = des.w1_mass_m;
    p.t1 = des.t1_m;
    p.E1 = des.E1_Pa;
    p.rho1 = des.rho1;
    p.delta = des.delta;
    p.P0 = des.P0_Pa;
    
    p.m1 = p.rho1 * p.L * p.w1_mass * p.t1;
    p.I1 = p.w1_stiff * p.t1^3 / 12;
    p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.K_couple = p.E1 * (p.w1_stiff * p.t1) / p.L;
    p.c1 = 2 * p.delta * sqrt(p.m1 * p.k1);
    
    p.S2 = des.S2_m2;
    p.rho2 = des.rho2;
    p.k2_base = des.k2_base_Npm;
    p.m2 = p.rho2 * p.S2 * p.L;
    p.c2 = 2 * p.delta * sqrt(p.m2 * p.k2_base);
    
    p.Gamma = (des.C_eff / 2) * 1000;
    
    % 平衡点計算 [file:1]
    coupling_eq = p.k1 / (2 * p.Gamma * p.K_couple);
    y_eq = p.Gamma * z_eq^2 + coupling_eq;
    y_natural = y_eq + (p.K_couple * coupling_eq) / p.k2_base;
    
    p.z_eq = z_eq; p.y_eq = y_eq; p.y_natural = y_natural;
    x0 = [z_eq; y_eq; 0; 0];
    
    % 3) 周波数スイープ
    [f_res, A_z] = find_resonance_with_amplitude(p, x0, cfg);
    
    out = struct(...
        'x_push_m', des.x_push_m, ...
        'w1_stiff_m', des.w1_stiff_m, ...
        'w1_mass_m', des.w1_mass_m, ...
        'L_m', des.L_m, ...
        'z_eq_mm', z_eq_mm, ...
        'f_res_Hz', f_res, ...
        'A_z_mm', A_z * 1000, ...
        'converged', true ...
    );
end

function [f_res, A_z_m] = find_resonance_with_amplitude(p, x0, cfg)
    % 粗スイープ
    amp_coarse = zeros(size(cfg.freq_coarse));
    for i = 1:length(cfg.freq_coarse)
        p.freq = cfg.freq_coarse(i);
        [~, x_sw] = ode45(@(t,x) eom_2dof(t,x,p,0), cfg.tspan_sweep, x0, cfg.ode_opt);
        z = x_sw(:,1);
        z_steady = z(round(end/2):end);
        amp_coarse(i) = (max(z_steady) - min(z_steady)) / 2;
    end
    [~, idx] = max(amp_coarse);
    f0 = cfg.freq_coarse(idx);
    
    % 細スイープ
    freq_fine = max(1, f0-cfg.freq_fine_half) : 1 : (f0+cfg.freq_fine_half);
    amp_fine = zeros(size(freq_fine));
    for i = 1:length(freq_fine)
        p.freq = freq_fine(i);
        [~, x_sw] = ode45(@(t,x) eom_2dof(t,x,p,0), cfg.tspan_sweep, x0, cfg.ode_opt);
        z = x_sw(:,1);
        z_steady = z(round(end/2):end);
        amp_fine(i) = (max(z_steady) - min(z_steady)) / 2;
    end
    [A_z_m, idx2] = max(amp_fine);
    f_res = freq_fine(idx2);
end

function dx = eom_2dof(t, x, p, t_sound_on)
    z = x(1); y = x(2); dz = x(3); dy = x(4);
    
    coupling_term = y - (p.Gamma * z^2);
    
    F_restore_z = p.k1 * z - p.K_couple * coupling_term * (2 * p.Gamma * z);
    F_restore_y = p.k2_base * (y - p.y_natural) + p.K_couple * coupling_term;
    
    if t >= t_sound_on
        ramp = min(1.0, (t - t_sound_on) / 0.1);
        F_sound = ramp * p.P0 * (p.L * p.w1_mass) * cos(2*pi*p.freq*t);
    else
        F_sound = 0;
    end
    
    ddz = (F_sound - p.c1*dz - F_restore_z) / p.m1;
    ddy = (0       - p.c2*dy - F_restore_y) / p.m2;
    
    dx = [dz; dy; ddz; ddy];
end

function plot_zeq_only_results(results)
    converged = [results.converged];
    res_ok = results(converged);
    
    if isempty(res_ok)
        warning('収束した設計点がありません');
        return;
    end
    
    z_eq_all = [res_ok.z_eq_mm];
    A_z_all = [res_ok.A_z_mm];
    f_res_all = [res_ok.f_res_Hz];
    
    [z_sorted, sort_idx] = sort(z_eq_all, 'descend');
    top20 = res_ok(sort_idx(1:min(20, length(res_ok))));
    
    figure('Color','w','Position',[100 100 1400 900]);
    
    subplot(2,3,1);
    bar([top20.z_eq_mm]);
    xlabel('順位'); ylabel('初期たわみ z_{eq} [mm]');
    title('トップ20設計（初期たわみ順）');
    grid on;
    
    subplot(2,3,2);
    histogram(z_eq_all, 30);
    xlabel('初期たわみ z_{eq} [mm]'); ylabel('設計点数');
    title('初期たわみの分布');
    grid on;
    
    subplot(2,3,3);
    scatter(z_eq_all, A_z_all, 30, f_res_all, 'filled', 'MarkerFaceAlpha', 0.6);
    colorbar; xlabel('初期たわみ z_{eq} [mm]'); ylabel('共振振幅 A_z [mm]');
    title('初期たわみ vs 共振振幅（色=共振周波数）');
    grid on;
    
    subplot(2,3,[4 5 6]);
    str_list = cell(min(20, length(top20)), 1);
    for i = 1:length(str_list)
        str_list{i} = sprintf('%2d: x_p=%.3f w_s=%.2f w_m=%5.1f L=%4.1f → z_eq=%.3f (A_z=%.3f, f=%5.1f Hz)', ...
            i, top20(i).x_push_m*1000, top20(i).w1_stiff_m*1000, ...
            top20(i).w1_mass_m*1000, top20(i).L_m*1000, ...
            top20(i).z_eq_mm, top20(i).A_z_mm, top20(i).f_res_Hz);
    end
    text(0.02, 0.98, str_list, 'FontSize', 9, 'VerticalAlignment', 'top', ...
         'Interpreter', 'none', 'FontName', 'FixedWidth');
    axis off;
    title('トップ20の詳細（単位: x_p,w_s,z_eq[mm], w_m,L[mm], f[Hz]）', 'FontSize', 11);
    
    fprintf('\n========== 初期たわみトップ20 ==========\n');
    fprintf('順位 | x_push  | w_stiff | w_mass | L    | z_eq   | A_z    | f_res \n');
    fprintf('-----|---------|---------|--------|------|--------|--------|-------\n');
    for i = 1:min(20, length(top20))
        fprintf('%4d | %6.3f  | %6.2f  | %5.1f  | %4.1f | %6.3f | %6.3f | %6.1f\n', ...
            i, top20(i).x_push_m*1000, top20(i).w1_stiff_m*1000, ...
            top20(i).w1_mass_m*1000, top20(i).L_m*1000, ...
            top20(i).z_eq_mm, top20(i).A_z_mm, top20(i).f_res_Hz);
    end
    fprintf('=========================================\n\n');
    
    best = top20(1);
    fprintf('【最大初期たわみ設計】\n');
    fprintf('  x_push   = %.3f mm\n', best.x_push_m*1000);
    fprintf('  w1_stiff = %.2f mm\n', best.w1_stiff_m*1000);
    fprintf('  w1_mass  = %.1f mm\n', best.w1_mass_m*1000);
    fprintf('  L        = %.1f mm\n', best.L_m*1000);
    fprintf('  → z_eq   = %.3f mm\n', best.z_eq_mm);
    fprintf('  （参考: A_z = %.3f mm, f_res = %.1f Hz）\n\n', best.A_z_mm, best.f_res_Hz);
end
