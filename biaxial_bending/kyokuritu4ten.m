% =========================================================
% 根本から6点の3次フィットのみで曲率を計算・プロット
% 片持はり / 部分固定
% 発表用：単位は mm，プロットは6点まで表示
% =========================================================

clear; clc; close all;

%% --- 1. データ入力（mmのまま使う） ---
x = [0, 0.41667, 0.83333, 1.25, 1.6667, 2.0833, 2.5, 2.9167, 3.3333, 3.75];

% 片持はり
y_cantilever = [ ...
     6.5564e-020 ...
     1.348e-006 ...
    -9.6879e-005 ...
    -3.544e-004 ...
    -7.8754e-004 ...
    -1.4003e-003 ...
    -2.1925e-003 ...
    -3.1608e-003 ...
    -4.3005e-003 ...
    -5.6065e-003];

% 部分固定
y_partial = [ ...
     4.4301e-019 ...
    -4.1087e-014 ...
    -1.2041e-003 ...
    -4.5133e-003 ...
    -8.9618e-003 ...
    -1.4106e-002 ...
    -1.9817e-002 ...
    -2.598e-002 ...
    -3.2521e-002 ...
    -3.9385e-002];

%% --- 2. 使用点数 ---
n_fit = 6;   % 根本から6点

%% --- 3. 計算 ---
result1 = calc_root_curvature_cubic(x, y_cantilever, n_fit, '片持はり');
result2 = calc_root_curvature_cubic(x, y_partial,    n_fit, '部分固定');

%% --- 4. 結果表示 ---
disp(' ');
disp('========== 根本から6点の3次フィットによる曲率 ==========');
print_result(result1);
print_result(result2);

%% --- 5. プロット ---
figure('Color','w','Name','6-point cubic fitting only (presentation, mm)');
tiledlayout(1,2);

plot_case_cubic(x, y_cantilever, n_fit, result1, '片持はり');
plot_case_cubic(x, y_partial,    n_fit, result2, '部分固定');

% =========================================================
% ローカル関数
% =========================================================

function result = calc_root_curvature_cubic(x, y, n_fit, name)

    % 3次フィット: y = ax^3 + bx^2 + cx + d
    p = polyfit(x(1:n_fit), y(1:n_fit), 3);

    b = p(2);
    c = p(3);

    yp0  = c;                              % y'(0)
    ypp0 = 2 * b;                          % y''(0)
    kappa0 = ypp0 / (1 + yp0^2)^(3/2);     % 根本曲率 [1/mm]
    R0 = 1 / abs(kappa0);                  % 曲率半径 [mm]

    result.name = name;
    result.p = p;
    result.yp0 = yp0;
    result.ypp0 = ypp0;
    result.kappa0 = kappa0;
    result.R0 = R0;
end

function print_result(r)
    fprintf('\n--- %s ---\n', r.name);
    fprintf('3次フィット y''''(0)   = %.6e [1/mm]\n', r.ypp0);
    fprintf('3次フィット 曲率      = %.6e [1/mm]\n', r.kappa0);
    fprintf('3次フィット 曲率半径  = %.6f [mm]\n', r.R0);
end

function plot_case_cubic(x, y, n_fit, r, fig_title)

    nexttile;
    hold on; grid on; box on;

    % 発表用：6点までのみ表示
    x_plot = x(1:n_fit);
    y_plot = y(1:n_fit);

    % 使用6点（凡例には出さない）
    plot(x_plot, y_plot, 'ro', ...
        'MarkerFaceColor', 'r', ...
        'MarkerSize', 6, ...
        'HandleVisibility', 'off');

    % 3次フィット曲線
    xfine = linspace(x_plot(1), x_plot(end), 200);
    yfit = polyval(r.p, xfine);
    plot(xfine, yfit, 'b-', ...
        'LineWidth', 2.0, ...
        'DisplayName', '6点3次フィット');

    title(sprintf('%s\n\\kappa(0)=%.6e [1/mm]', fig_title, r.kappa0), 'FontSize', 12);
    xlabel('x [mm]');
    ylabel('y [mm]');

    % 軸範囲を固定
    xlim([0 2.2]);
    ylim([-1.45e-2 1e-3]);

    legend('Location', 'best');
end