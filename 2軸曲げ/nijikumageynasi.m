% ==========================================
% y方向の境界条件を外して、x方向たわみをy方向に一様展開
% ==========================================
clear; clc; close all;

% 1. 材料・形状パラメータの定義
E = 3.45e9;      % ヤング率 [Pa]
nu = 0.39;       % ポアソン比 [-]
h = 0.1e-3;      % 厚さ [m]
lx = 20e-3;      % x方向の長さ [m]
ly = 18e-3;      % y方向の長さ [m]
q = 2;           % 等分布荷重 [Pa]

% 2. 係数の計算
n = 1;
alpha = (n * pi) / ly;
a = alpha * lx;
qn = (4 * q) / (n * pi);
D_flex = (E * h^3) / (12 * (1 - nu^2));
K = qn / (D_flex * alpha^4);

% 3. 係数行列 M と右辺ベクトル b
M = zeros(4, 4);

% W(0) = 0
M(1, :) = [1, 0, 0, 0];

% W'(0) = 0
M(2, :) = [0, 1, 1, 0];

% M_x(lx) = 0
M(3, :) = [(1-nu)*cosh(a), ...
           (1-nu)*sinh(a), ...
           2*sinh(a)+a*(1-nu)*cosh(a), ...
           2*cosh(a)+a*(1-nu)*sinh(a)];

% V_x(lx) = 0
M(4, :) = [(nu-1)*sinh(a), ...
           (nu-1)*cosh(a), ...
           (1+nu)*cosh(a)+a*(nu-1)*sinh(a), ...
           (1+nu)*sinh(a)+a*(nu-1)*cosh(a)];

b = [-K; 0; nu*K; 0];

% 4. 連立方程式を解いて係数 A, B, C, D を求める
coeffs = M \ b;
A = coeffs(1);
B = coeffs(2);
C = coeffs(3);
D = coeffs(4);

% 5. 空間メッシュの作成
x_vec = linspace(0, lx, 100);
y_vec = linspace(0, ly, 100);
[X, Y] = meshgrid(x_vec, y_vec);

% 6. 静的たわみ w(x,y) の計算
Z = alpha * x_vec;
W_x = A*cosh(Z) + B*sinh(Z) + C.*Z.*cosh(Z) + D.*Z.*sinh(Z) + K;

% y方向には依存させず、一様に展開
w = repmat(W_x, length(y_vec), 1);

% ==========================================
% たわみだけをプロット
% ==========================================
figure('Position', [100, 100, 700, 500]);
surf(X*1000, Y*1000, w*1000, 'EdgeColor', 'none');
colormap(parula);
colorbar;

title('静的たわみ分布 w [mm]');
xlabel('x (固定端→自由端) [mm]');
ylabel('y [mm]');
zlabel('たわみ [mm]');
view(-45, 30);