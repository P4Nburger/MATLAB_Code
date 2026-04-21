clear; clc; close all;

% ===== 入力値（単位:mm）=====
x1 = 0.0;   y1 = 0.0;
x2 = 1;  y2 = -2.364e-3;
x3 = 2; y3 = -1.3026e-2;


% ===== 3点から曲率半径・曲率を計算 =====
A2 = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);

s1 = sqrt((x2-x3)^2 + (y2-y3)^2);
s2 = sqrt((x1-x3)^2 + (y1-y3)^2);
s3 = sqrt((x1-x2)^2 + (y1-y2)^2);

eps_val = 1e-12;

if abs(A2) < eps_val
    R = Inf;          % 曲率半径 [mm]
    kappa = 0;        % 曲率 [1/mm]
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

% ===== 結果表示 =====
fprintf('R = %.6f mm\n', R);
fprintf('kappa = %.6f 1/mm\n', kappa);

% ===== プロット =====
figure;
hold on; grid on; axis equal;

% 3点をプロット
plot([x1 x2 x3], [y1 y2 y3], 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot([x1 x2 x3], [y1 y2 y3], 'k--', 'LineWidth', 1);

% 円をプロット（一直線でない場合）
if isfinite(R)
    th = linspace(0, 2*pi, 400);
    xcirc = xc + R*cos(th);
    ycirc = yc + R*sin(th);
    plot(xcirc, ycirc, 'b-', 'LineWidth', 1.5);
    plot(xc, yc, 'bx', 'MarkerSize', 10, 'LineWidth', 2);
    title(sprintf('R = %.4f mm,  kappa = %.6f 1/mm', R, kappa));
    legend('3 points', 'Polyline', 'Fitted circle', 'Circle center', 'Location', 'best');
else
    title('3 points are collinear: R = Inf, kappa = 0');
    legend('3 points', 'Polyline', 'Location', 'best');
end

xlabel('x [mm]');
ylabel('y [mm]');