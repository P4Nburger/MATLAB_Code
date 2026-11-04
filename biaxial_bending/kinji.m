% =========================================================
% Kirchhoff-Love plate bending by Ritz method
% No separation of variables
% x=0: clamped  -> w=0, dw/dx=0
% Other edges: natural boundaries (free unless additional terms are imposed)
% =========================================================
clear; clc; close all;

%% 1) Material / geometry / load
E  = 3.45e9;      % Young's modulus [Pa]
nu = 0.39;        % Poisson ratio [-]
h  = 0.1e-3;      % thickness [m]  <-- adjust if needed
lx = 20e-3;       % length in x [m]
ly = 18e-3;       % length in y [m]
q  = 2;           % uniform transverse load [Pa]

D = E*h^3/(12*(1-nu^2));   % flexural rigidity

%% 2) Ritz basis
% Dimensionless coordinates:
% xi  = x/lx        in [0,1]
% eta = 2*y/ly - 1  in [-1,1]
%
% Symmetric basis in y for uniform symmetric loading.
% All terms include xi^2, so at x=0:
%   w = 0
%   dw/dx = 0
%
% phi_k = xi^px * eta^py
basis_exp = [
    2 0
    3 0
    4 0
    2 2
    3 2
    4 2
    2 4
    3 4
];

nb = size(basis_exp,1);

%% 3) Build stiffness matrix K and load vector f
K = zeros(nb, nb);
f = zeros(nb, 1);

for i = 1:nb
    pi_x = basis_exp(i,1);
    pi_y = basis_exp(i,2);

    % load vector
    f(i) = integral2(@(x,y) q .* phi_xy(x,y,pi_x,pi_y,lx,ly), ...
                     0, lx, 0, ly, 'AbsTol',1e-10, 'RelTol',1e-8);

    for j = 1:nb
        pj_x = basis_exp(j,1);
        pj_y = basis_exp(j,2);

        integrand = @(x,y) ...
            d2phidx2(x,y,pj_x,pj_y,lx,ly).*d2phidx2(x,y,pi_x,pi_y,lx,ly) + ...
            d2phidy2(x,y,pj_x,pj_y,lx,ly).*d2phidy2(x,y,pi_x,pi_y,lx,ly) + ...
            nu .* d2phidx2(x,y,pj_x,pj_y,lx,ly).*d2phidy2(x,y,pi_x,pi_y,lx,ly) + ...
            nu .* d2phidy2(x,y,pj_x,pj_y,lx,ly).*d2phidx2(x,y,pi_x,pi_y,lx,ly) + ...
            2*(1-nu).*d2phidxdy(x,y,pj_x,pj_y,lx,ly).*d2phidxdy(x,y,pi_x,pi_y,lx,ly);

        K(i,j) = D * integral2(integrand, 0, lx, 0, ly, ...
                               'AbsTol',1e-10, 'RelTol',1e-8);
    end
end

%% 4) Solve Ritz coefficients
a = K \ f;

disp('Ritz coefficients a = ');
disp(a);

%% 5) Evaluate displacement and moments on grid
nx = 101;
ny = 101;
x_vec = linspace(0, lx, nx);
y_vec = linspace(0, ly, ny);
[X, Y] = meshgrid(x_vec, y_vec);

w   = zeros(size(X));
wxx = zeros(size(X));
wyy = zeros(size(X));
wxy = zeros(size(X));

for k = 1:nb
    px = basis_exp(k,1);
    py = basis_exp(k,2);

    w   = w   + a(k)*phi_xy(X,Y,px,py,lx,ly);
    wxx = wxx + a(k)*d2phidx2(X,Y,px,py,lx,ly);
    wyy = wyy + a(k)*d2phidy2(X,Y,px,py,lx,ly);
    wxy = wxy + a(k)*d2phidxdy(X,Y,px,py,lx,ly);
end

% Moment-curvature relations
Mx  = -D*(wxx + nu*wyy);
My  = -D*(wyy + nu*wxx);
Mxy = -D*(1-nu)*wxy;

%% 6) Report max values
fprintf('max deflection = %.6e m\n', max(abs(w(:))));
fprintf('max |Mx|       = %.6e N\n', max(abs(Mx(:))));
fprintf('max |My|       = %.6e N\n', max(abs(My(:))));
fprintf('max |Mxy|      = %.6e N\n', max(abs(Mxy(:))));

%% 7) Plot deflection
figure('Position',[100 100 760 560]);
surf(X*1e3, Y*1e3, w*1e3, 'EdgeColor','none');
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('w [mm]');
title('Ritz法によるたわみ分布');
colormap(parula); colorbar; view(-40,30);

%% 8) Plot moments
figure('Position',[120 120 1200 360]);

subplot(1,3,1);
surf(X*1e3, Y*1e3, Mx, 'EdgeColor','none');
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('Mx [N]');
title('M_x'); colorbar; view(-40,30);

subplot(1,3,2);
surf(X*1e3, Y*1e3, My, 'EdgeColor','none');
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('My [N]');
title('M_y'); colorbar; view(-40,30);

subplot(1,3,3);
surf(X*1e3, Y*1e3, Mxy, 'EdgeColor','none');
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('M_{xy} [N]');
title('M_{xy}'); colorbar; view(-40,30);

% =========================================================
% Local functions
% =========================================================
function val = phi_xy(x,y,px,py,lx,ly)
    xi  = x ./ lx;
    eta = 2*y ./ ly - 1;
    val = xi.^px .* eta.^py;
end

function val = d2phidx2(x,y,px,py,lx,ly)
    xi  = x ./ lx;
    eta = 2*y ./ ly - 1;
    if px < 2
        val = zeros(size(x));
    else
        val = (px*(px-1)/lx^2) .* xi.^(px-2) .* eta.^py;
    end
end

function val = d2phidy2(x,y,px,py,lx,ly)
    xi  = x ./ lx;
    eta = 2*y ./ ly - 1;
    if py < 2
        val = zeros(size(x));
    else
        val = (py*(py-1)*(2/ly)^2) .* xi.^px .* eta.^(py-2);
    end
end

function val = d2phidxdy(x,y,px,py,lx,ly)
    xi  = x ./ lx;
    eta = 2*y ./ ly - 1;
    if px < 1 || py < 1
        val = zeros(size(x));
    else
        val = (px*py*(1/lx)*(2/ly)) .* xi.^(px-1) .* eta.^(py-1);
    end
end