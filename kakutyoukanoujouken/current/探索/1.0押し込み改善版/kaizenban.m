clear; clc; close all;

delete(gcp('nocreate'));
pool = parpool('Threads', 20);

fprintf('=======================================================\n');
fprintf('   Extended Parameter Search: x_push up to 1.0 mm\n');
fprintf('   Enhanced Analysis with Top 20 Verification\n');
fprintf('=======================================================\n\n');

% Extended search space
ss = struct(...
    'x_push_m', linspace(0.02e-3, 1.0e-3, 20), ...
    'w1_stiff_m', [0.1 0.5 1 2 5]*1e-3, ...
    'w1_mass_m', [10 20 40 60 80 100]*1e-3, ...
    'L_m', [15 18 21 25]*1e-3 ...
);

fx = struct(...
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

cfg = struct(...
    'freq_coarse', 10:10:200, ...
    'freq_fine_half', 10, ...
    'tspan_sweep', [0 2.0], ...
    'ode_opt', odeset('RelTol',1e-6,'AbsTol',1e-9) ...
);

n_tot = prod(structfun(@length, ss));
fprintf('Total design points: %d (20x5x6x4)\n', n_tot);
fprintf('Parallel workers: %d\n', pool.NumWorkers);
fprintf('Estimated time: %.1f min\n', 6*2400/960);
fprintf('Start time: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

% Create output directory
out_dir = sprintf('results_%s', datestr(now, 'yyyymmdd_HHMMSS'));
mkdir(out_dir);
fprintf('Output directory: %s\n\n', out_dir);

tic;
res = run_par_with_progress(ss, fx, cfg, n_tot);
t_elap = toc;

% Save results
mat_file = fullfile(out_dir, 'zeq_extended_search.mat');
save(mat_file, 'res', 'ss', 'fx', 'cfg', 't_elap');
fprintf('\n[SAVED] Main results: %s\n', mat_file);

% Plot and analyze with enhanced visualization
[top_designs, fig_handle] = plot_res_enhanced(res, out_dir);

% Export top 20 to CSV
csv_file = fullfile(out_dir, 'top20_designs.csv');
export_top_designs(top_designs, csv_file);

% Automatic snapthrough verification for TOP 20
n_verify = min(20, length(top_designs));
fprintf('\n=== Snapthrough Verification (Top %d) ===\n', n_verify);
verify_snapthrough(top_designs(1:n_verify), fx, cfg, out_dir);

fprintf('\n=======================================================\n');
fprintf('   ANALYSIS COMPLETED\n');
fprintf('=======================================================\n');
fprintf('Total time: %.1f min\n', t_elap/60);
fprintf('End time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('All results saved in: %s\n', out_dir);
fprintf('=======================================================\n');

%% ========== FUNCTIONS ==========

function res = run_par_with_progress(ss, fx, cfg, n_tot)
    [Xp, Ws, Wm, Ls] = ndgrid(ss.x_push_m, ss.w1_stiff_m, ss.w1_mass_m, ss.L_m);
    Xp = Xp(:); Ws = Ws(:); Wm = Wm(:); Ls = Ls(:);
    
    z_arr = zeros(n_tot, 1);
    f_arr = zeros(n_tot, 1);
    a_arr = zeros(n_tot, 1);
    c_arr = false(n_tot, 1);
    
    fprintf('Computing %d designs...\n', n_tot);
    t_start = tic;
    
    parfor i = 1:n_tot
        d = struct('x_push_m', Xp(i), 'w1_stiff_m', Ws(i), ...
                   'w1_mass_m', Wm(i), 'L_m', Ls(i), ...
                   't1_m', fx.t1_m, 'E1_Pa', fx.E1_Pa, 'rho1', fx.rho1, ...
                   'P0_Pa', fx.P0_Pa, 'delta', fx.delta, 'C_eff', fx.C_eff, ...
                   'S2_m2', fx.S2_m2, 'rho2', fx.rho2, 'k2_base_Npm', fx.k2_base_Npm);
        try
            o = eval_d(d, cfg);
            z_arr(i) = o.z_eq_mm;
            f_arr(i) = o.f_res_Hz;
            a_arr(i) = o.A_z_mm;
            c_arr(i) = true;
        catch ME
            z_arr(i) = NaN;
            f_arr(i) = NaN;
            a_arr(i) = NaN;
            c_arr(i) = false;
        end
    end
    
    elapsed = toc(t_start);
    fprintf('Computation completed in %.1f min\n', elapsed/60);
    
    % Build result structure
    res = struct('x_push_m', {}, 'w1_stiff_m', {}, 'w1_mass_m', {}, 'L_m', {}, ...
                 'z_eq_mm', {}, 'f_res_Hz', {}, 'A_z_mm', {}, 'converged', {});
    for i = 1:n_tot
        res(i).x_push_m = Xp(i);
        res(i).w1_stiff_m = Ws(i);
        res(i).w1_mass_m = Wm(i);
        res(i).L_m = Ls(i);
        res(i).z_eq_mm = z_arr(i);
        res(i).f_res_Hz = f_arr(i);
        res(i).A_z_mm = a_arr(i);
        res(i).converged = c_arr(i);
    end
    
    fprintf('Converged: %d / %d (%.1f%%)\n', sum(c_arr), n_tot, 100*sum(c_arr)/n_tot);
end

function o = eval_d(d, cfg)
    z_mm = sqrt(2 * (d.x_push_m*1000) / d.C_eff);
    z = z_mm / 1000;
    p.L = d.L_m; p.w1_stiff = d.w1_stiff_m; p.w1_mass = d.w1_mass_m;
    p.t1 = d.t1_m; p.E1 = d.E1_Pa; p.rho1 = d.rho1;
    p.delta = d.delta; p.P0 = d.P0_Pa;
    p.m1 = p.rho1 * p.L * p.w1_mass * p.t1;
    p.I1 = p.w1_stiff * p.t1^3 / 12;
    p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.Kc = p.E1 * (p.w1_stiff * p.t1) / p.L;
    p.c1 = 2 * p.delta * sqrt(p.m1 * p.k1);
    p.S2 = d.S2_m2; p.rho2 = d.rho2; p.k2 = d.k2_base_Npm;
    p.m2 = p.rho2 * p.S2 * p.L;
    p.c2 = 2 * p.delta * sqrt(p.m2 * p.k2);
    p.Gm = (d.C_eff / 2) * 1000;
    ceq = p.k1 / (2 * p.Gm * p.Kc);
    yeq = p.Gm * z^2 + ceq;
    ynat = yeq + (p.Kc * ceq) / p.k2;
    p.z_eq = z; p.y_eq = yeq; p.y_nat = ynat;
    x0 = [z; yeq; 0; 0];
    [fr, Az] = find_r(p, x0, cfg);
    o = struct('z_eq_mm', z_mm, 'f_res_Hz', fr, 'A_z_mm', Az*1000);
end

function [fr, Az] = find_r(p, x0, cfg)
    ac = zeros(size(cfg.freq_coarse));
    for i = 1:length(cfg.freq_coarse)
        p.freq = cfg.freq_coarse(i);
        [~, xs] = ode45(@(t,x) eom(t,x,p,0), cfg.tspan_sweep, x0, cfg.ode_opt);
        z = xs(:,1); zs = z(round(end/2):end);
        ac(i) = (max(zs) - min(zs)) / 2;
    end
    [~, ix] = max(ac); f0 = cfg.freq_coarse(ix);
    ff = max(1, f0-cfg.freq_fine_half) : 1 : (f0+cfg.freq_fine_half);
    af = zeros(size(ff));
    for i = 1:length(ff)
        p.freq = ff(i);
        [~, xs] = ode45(@(t,x) eom(t,x,p,0), cfg.tspan_sweep, x0, cfg.ode_opt);
        z = xs(:,1); zs = z(round(end/2):end);
        af(i) = (max(zs) - min(zs)) / 2;
    end
    [Az, ix2] = max(af); fr = ff(ix2);
end

function dx = eom(t, x, p, t0)
    z = x(1); y = x(2); dz = x(3); dy = x(4);
    ct = y - (p.Gm * z^2);
    Fz = p.k1 * z - p.Kc * ct * (2 * p.Gm * z);
    Fy = p.k2 * (y - p.y_nat) + p.Kc * ct;
    if t >= t0
        rp = min(1.0, (t - t0) / 0.1);
        Fs = rp * p.P0 * (p.L * p.w1_mass) * cos(2*pi*p.freq*t);
    else
        Fs = 0;
    end
    ddz = (Fs - p.c1*dz - Fz) / p.m1;
    ddy = (0  - p.c2*dy - Fy) / p.m2;
    dx = [dz; dy; ddz; ddy];
end

function [top_designs, fig_handle] = plot_res_enhanced(res, out_dir)
    cv = [res.converged]; r = res(cv);
    if isempty(r), fprintf('No converged results\n'); top_designs = []; fig_handle = []; return; end
    
    za = [r.z_eq_mm]; aa = [r.A_z_mm]; fa = [r.f_res_Hz]; 
    xp = [r.x_push_m]*1000; ws = [r.w1_stiff_m]*1000;
    wm = [r.w1_mass_m]*1000; Ls = [r.L_m]*1000;
    ratio_all = aa ./ za;
    
    [~, ix] = sort(aa, 'descend');
    top_designs = r(ix(1:min(20, length(r))));
    ratio_top = [top_designs.A_z_mm] ./ [top_designs.z_eq_mm];
    
    % ========== Enhanced figure: 8 subplots ==========
    fig_handle = figure('Color','w','Position',[50 50 1800 1000]);
    
    % (1) Top 20 with A_z/z_eq ratio
    subplot(2,4,1);
    bar(ratio_top); hold on;
    yline(1.0, 'r--', 'LineWidth', 2);
    xlabel('Rank'); ylabel('A_z / z_{eq}');
    title('Top 20: Snapthrough Indicator'); 
    legend('Ratio', 'Threshold=1.0', 'Location','northeast');
    grid on;
    
    % (2) Distribution of A_z/z_eq ratio
    subplot(2,4,2);
    histogram(ratio_all, 40); hold on;
    xline(1.0, 'r--', 'LineWidth', 2);
    xlabel('A_z / z_{eq}'); ylabel('Count');
    title('Distribution (snapthrough if >1.0)');
    n_snap = sum(ratio_all > 1.0);
    n_total = length(ratio_all);
    text(1.5, max(ylim)*0.8, sprintf('%.1f%% > 1.0', 100*n_snap/n_total), ...
        'FontSize', 10, 'FontWeight', 'bold');
    grid on;
    
    % (3) Design space map
    subplot(2,4,3);
    scatter(xp, za, 30, ratio_all, 'filled', 'MarkerFaceAlpha', 0.7);
    colorbar; caxis([0 2]); colormap(gca, jet);
    hold on; 
    x_theory = linspace(0, 1, 100);
    z_theory = sqrt(2*x_theory/0.2794);
    plot(x_theory, z_theory, 'k--', 'LineWidth', 1.5);
    xlabel('x_{push} [mm]'); ylabel('z_{eq} [mm]');
    title('Design Space (color=A_z/z_{eq})');
    legend('Data', 'Theory z_{eq}', 'Location','northwest');
    grid on;
    
    % (4) x_push vs A_z
    subplot(2,4,4);
    scatter(xp, aa, 30, ratio_all, 'filled', 'MarkerFaceAlpha', 0.6);
    colorbar; caxis([0 2]);
    xlabel('x_{push} [mm]'); ylabel('A_z [mm]');
    title('x_{push} effect (color=ratio)');
    grid on;
    
    % (5) w1_stiff vs A_z
    subplot(2,4,5);
    scatter(ws, aa, 30, ratio_all, 'filled', 'MarkerFaceAlpha', 0.6);
    set(gca, 'XScale', 'log');
    colorbar; caxis([0 2]);
    xlabel('w1_{stiff} [mm]'); ylabel('A_z [mm]');
    title('w1_{stiff} effect (log scale)');
    grid on;
    
    % (6) w1_mass vs A_z
    subplot(2,4,6);
    scatter(wm, aa, 30, ratio_all, 'filled', 'MarkerFaceAlpha', 0.6);
    colorbar; caxis([0 2]);
    xlabel('w1_{mass} [mm]'); ylabel('A_z [mm]');
    title('w1_{mass} effect');
    grid on;
    
    % (7) L vs A_z
    subplot(2,4,7);
    scatter(Ls, aa, 30, ratio_all, 'filled', 'MarkerFaceAlpha', 0.6);
    colorbar; caxis([0 2]);
    xlabel('L [mm]'); ylabel('A_z [mm]');
    title('L (span) effect');
    grid on;
    
    % (8) f_res vs A_z
    subplot(2,4,8);
    scatter(fa, aa, 30, ratio_all, 'filled', 'MarkerFaceAlpha', 0.6);
    colorbar; caxis([0 2]);
    xlabel('f_{res} [Hz]'); ylabel('A_z [mm]');
    title('Resonance frequency effect');
    grid on;
    
    % Save
    fig_file = fullfile(out_dir, 'results_summary_enhanced.png');
    saveas(fig_handle, fig_file);
    fprintf('[SAVED] Enhanced figure: %s\n', fig_file);
    
    % Print results
    fprintf('\n========== TOP 20 (Sorted by A_z) ==========\n');
    fprintf('Rank | x_push | w_stiff | w_mass | L    | A_z    | z_eq   | f_res | Az/zeq\n');
    fprintf('-----|--------|---------|--------|------|--------|--------|-------|-------\n');
    for i = 1:min(20, length(top_designs))
        fprintf('%4d | %6.3f | %7.2f | %6.1f | %4.1f | %6.3f | %6.3f | %6.1f | %.3f\n', ...
            i, top_designs(i).x_push_m*1000, top_designs(i).w1_stiff_m*1000, ...
            top_designs(i).w1_mass_m*1000, top_designs(i).L_m*1000, ...
            top_designs(i).A_z_mm, top_designs(i).z_eq_mm, ...
            top_designs(i).f_res_Hz, ratio_top(i));
    end
    fprintf('============================================\n\n');
    
    % Best design summary
    best = top_designs(1);
    fprintf('BEST DESIGN (Rank 1):\n');
    fprintf('  x_push   = %.3f mm\n', best.x_push_m*1000);
    fprintf('  w1_stiff = %.2f mm\n', best.w1_stiff_m*1000);
    fprintf('  w1_mass  = %.1f mm\n', best.w1_mass_m*1000);
    fprintf('  L        = %.1f mm\n', best.L_m*1000);
    fprintf('  A_z      = %.3f mm  <-- MAX RESONANCE AMPLITUDE\n', best.A_z_mm);
    fprintf('  z_eq     = %.3f mm\n', best.z_eq_mm);
    fprintf('  f_res    = %.1f Hz\n', best.f_res_Hz);
    fprintf('  A_z/z_eq = %.3f  ', ratio_top(1));
    if ratio_top(1) > 1.0
        fprintf('<-- SNAPTHROUGH CONFIRMED\n\n');
    else
        fprintf('<-- No snapthrough\n\n');
    end
    
    fprintf('Statistics:\n');
    fprintf('  Designs with A_z/z_eq > 1.0: %d / %d (%.1f%%)\n', ...
        n_snap, n_total, 100*n_snap/n_total);
    fprintf('  Max A_z/z_eq ratio: %.3f\n', max(ratio_all));
    fprintf('  Mean A_z/z_eq ratio: %.3f\n\n', mean(ratio_all));
end

function export_top_designs(top_designs, csv_file)
    fid = fopen(csv_file, 'w');
    fprintf(fid, 'Rank,x_push_mm,w1_stiff_mm,w1_mass_mm,L_mm,A_z_mm,z_eq_mm,f_res_Hz,Az_zeq_ratio\n');
    for i = 1:length(top_designs)
        ratio = top_designs(i).A_z_mm / top_designs(i).z_eq_mm;
        fprintf(fid, '%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%.4f\n', ...
            i, top_designs(i).x_push_m*1000, top_designs(i).w1_stiff_m*1000, ...
            top_designs(i).w1_mass_m*1000, top_designs(i).L_m*1000, ...
            top_designs(i).A_z_mm, top_designs(i).z_eq_mm, ...
            top_designs(i).f_res_Hz, ratio);
    end
    fclose(fid);
    fprintf('[SAVED] CSV: %s\n', csv_file);
end

function verify_snapthrough(top_designs, fx, cfg, out_dir)
    n_designs = length(top_designs);
    snap_count = 0;
    
    for i = 1:n_designs
        d = top_designs(i);
        fprintf('\n--- Rank %d/%d: xp=%.3f ws=%.2f wm=%.1f L=%.1f ---\n', ...
            i, n_designs, d.x_push_m*1000, d.w1_stiff_m*1000, ...
            d.w1_mass_m*1000, d.L_m*1000);
        
        % Build parameters
        z_mm = d.z_eq_mm; z = z_mm / 1000;
        p.L = d.L_m; p.w1_stiff = d.w1_stiff_m; p.w1_mass = d.w1_mass_m;
        p.t1 = fx.t1_m; p.E1 = fx.E1_Pa; p.rho1 = fx.rho1;
        p.delta = fx.delta; p.P0 = fx.P0_Pa;
        p.m1 = p.rho1 * p.L * p.w1_mass * p.t1;
        p.I1 = p.w1_stiff * p.t1^3 / 12;
        p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
        p.Kc = p.E1 * (p.w1_stiff * p.t1) / p.L;
        p.c1 = 2 * p.delta * sqrt(p.m1 * p.k1);
        p.S2 = fx.S2_m2; p.rho2 = fx.rho2; p.k2 = fx.k2_base_Npm;
        p.m2 = p.rho2 * p.S2 * p.L;
        p.c2 = 2 * p.delta * sqrt(p.m2 * p.k2);
        p.Gm = (fx.C_eff / 2) * 1000;
        ceq = p.k1 / (2 * p.Gm * p.Kc);
        yeq = p.Gm * z^2 + ceq;
        ynat = yeq + (p.Kc * ceq) / p.k2;
        p.z_eq = z; p.y_eq = yeq; p.y_nat = ynat; p.freq = d.f_res_Hz;
        
        x0 = [z; yeq; 0; 0];
        
        % Time history (5 seconds)
        t_on = 0.5;
        [t, x] = ode45(@(t,x) eom(t,x,p,t_on), [0 5], x0, cfg.ode_opt);
        
        % Snapthrough detection
        z_disp = x(:,1);
        threshold = 0.05 * abs(z);  % 5% threshold
        crosses_zero = any(z_disp > threshold) && any(z_disp < -threshold);
        
        if crosses_zero
            fprintf('  >> SNAPTHROUGH CONFIRMED (z crosses zero)\n');
            snap_count = snap_count + 1;
        else
            fprintf('  >> No snapthrough (z stays on one side)\n');
        end
        
        fprintf('  z_eq=%.3f mm, max(z)=%.3f mm, min(z)=%.3f mm\n', ...
            z*1000, max(z_disp)*1000, min(z_disp)*1000);
        fprintf('  A_z=%.3f mm, A_z/z_eq=%.3f\n', d.A_z_mm, d.A_z_mm/d.z_eq_mm);
        
        % Save plot
        fig = figure('Visible','off','Position',[100 100 1200 500]);
        subplot(1,2,1);
        plot(t, z_disp*1000, 'b', 'LineWidth', 1.5); hold on;
        xline(t_on, 'k--', 'LineWidth', 1.5);
        yline(0, 'r--', 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('z [mm]');
        title(sprintf('Rank %d: z displacement (A_z/z_{eq}=%.3f)', i, d.A_z_mm/d.z_eq_mm));
        grid on;
        if crosses_zero
            text(2.5, max(z_disp)*1000*0.8, 'SNAPTHROUGH!', ...
                'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
        end
        
        subplot(1,2,2);
        plot(t, (x(:,2)-yeq)*1000, 'r', 'LineWidth', 1.5); hold on;
        xline(t_on, 'k--', 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Frame expansion [mm]');
        title(sprintf('Rank %d: Frame expansion', i));
        grid on;
        
        png_file = fullfile(out_dir, sprintf('snapthrough_rank%02d.png', i));
        saveas(fig, png_file);
        close(fig);
        fprintf('  [SAVED] %s\n', png_file);
    end
    
    fprintf('\n========== SNAPTHROUGH SUMMARY ==========\n');
    fprintf('Verified: %d designs\n', n_designs);
    fprintf('Snapthrough confirmed: %d / %d (%.1f%%)\n', ...
        snap_count, n_designs, 100*snap_count/n_designs);
    fprintf('=========================================\n');
end
