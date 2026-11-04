clear; clc; close all;

% ===== 入力値（単位:mm）=====
x1 = 0.0;   y1 = 0.0;
x2 = 1;  y2 = -2.3644e-3;
x3 = 2; y3 = -1.3026e-2;


% ===== 3点から曲率半径・曲率を計算 =====
A2 = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);

s1 = sqrt((x2-x3)^2 + (y2-y3)^2);
s2 = sqrt((x1-x3)^2 + (y1-y3)^2);
s3 = sqrt((x1-x2)^2 + (y1-y2)^2);

eps_val = 1e-12;

if abs(A2) < eps_val
    R = Inf;
    kappa = 0;
    xc = NaN;
    yc = NaN;
else
    A = abs(A2)/2;
    R = (s1*s2*s3)/(4*A);
    kappa = 1/R;

    D = 2*(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
    xc = ((x1^2+y1^2)*(y2-y3) + (x2^2+y2^2)*(y3-y1) + (x3^2+y3^2)*(y1-y2)) / D;
    yc = ((x1^2+y1^2)*(x3-x2) + (x2^2+y2^2)*(x1-x3) + (x3^2+y3^2)*(x2-x1)) / D;
end

fprintf('R = %.6f mm\n', R);
fprintf('kappa = %.6f 1/mm\n', kappa);

figure;
hold on; grid on; axis equal;

% 3点
plot([x1 x2 x3], [y1 y2 y3], 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

if isfinite(R)
    % ===== 円全体 =====
    th_circle = linspace(0, 2*pi, 2000);
    x_circle = xc + R*cos(th_circle);
    y_circle = yc + R*sin(th_circle);
    plot(x_circle, y_circle, 'c--', 'LineWidth', 1.2);

    % ===== 3点を通る円弧 =====
    t1 = atan2(y1 - yc, x1 - xc);
    t2 = atan2(y2 - yc, x2 - xc);
    t3 = atan2(y3 - yc, x3 - xc);

    if t1 < 0, t1 = t1 + 2*pi; end
    if t2 < 0, t2 = t2 + 2*pi; end
    if t3 < 0, t3 = t3 + 2*pi; end

    t3_ccw = t3;
    if t3_ccw < t1
        t3_ccw = t3_ccw + 2*pi;
    end

    t2_test = t2;
    if t2_test < t1
        t2_test = t2_test + 2*pi;
    end

    if t1 <= t2_test && t2_test <= t3_ccw
        th_arc = linspace(t1, t3_ccw, 1000);
    else
        if t1 < t3
            t1 = t1 + 2*pi;
        end
        th_arc = linspace(t1, t3, 1000);
    end

    x_arc = xc + R*cos(th_arc);
    y_arc = yc + R*sin(th_arc);
    plot(x_arc, y_arc, 'b-', 'LineWidth', 2.5);

    % 円の中心
    plot(xc, yc, 'bx', 'MarkerSize', 10, 'LineWidth', 2);

    title(sprintf('R = %.4f mm,  kappa = %.6f 1/mm', R, kappa));
    legend('3 points', 'Full circle', 'Fitted arc', 'Circle center', 'Location', 'best');
else
    title('3 points are collinear: R = Inf, kappa = 0');
    legend('3 points', 'Location', 'best');
end

xlabel('x [mm]');
ylabel('y [mm]');