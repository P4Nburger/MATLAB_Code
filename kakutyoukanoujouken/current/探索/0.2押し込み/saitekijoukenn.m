% Best design parameters (Rank 16)
des_best = struct(...
    'x_push_m', 0.200e-3, ...
    'w1_stiff_m', 0.10e-3, ...
    'w1_mass_m', 60.0e-3, ...
    'L_m', 15.0e-3, ...
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

% Parameters calculation
p = struct();
p.L = des_best.L_m;
p.w1_stiff = des_best.w1_stiff_m;
p.w1_mass = des_best.w1_mass_m;
p.t1 = des_best.t1_m;
p.E1 = des_best.E1_Pa;
p.rho1 = des_best.rho1;
p.delta = des_best.delta;
p.P0 = des_best.P0_Pa;

p.m1 = p.rho1 * p.L * p.w1_mass * p.t1;
p.I1 = p.w1_stiff * p.t1^3 / 12;
p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
p.Kc = p.E1 * (p.w1_stiff * p.t1) / p.L;
p.c1 = 2 * p.delta * sqrt(p.m1 * p.k1);

p.S2 = des_best.S2_m2;
p.rho2 = des_best.rho2;
p.k2 = des_best.k2_base_Npm;
p.m2 = p.rho2 * p.S2 * p.L;
p.c2 = 2 * p.delta * sqrt(p.m2 * p.k2);

p.Gm = (des_best.C_eff / 2) * 1000;

% Initial deflection
z_mm = sqrt(2 * (des_best.x_push_m*1000) / des_best.C_eff);
z_eq = z_mm / 1000;

% Equilibrium calculation
ceq = p.k1 / (2 * p.Gm * p.Kc);
y_eq = p.Gm * z_eq^2 + ceq;
y_nat = y_eq + (p.Kc * ceq) / p.k2;

p.z_eq = z_eq;
p.y_eq = y_eq;
p.y_nat = y_nat;
p.freq = 33;  % Resonance frequency

% Initial conditions
x0 = [z_eq; y_eq; 0; 0];

% Time simulation
t_sound_on = 0.5;  % Sound starts at 0.5s
tspan = [0 5];  % 5 seconds
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);

[t, x] = ode45(@(t,x) eom_expansion(t,x,p,t_sound_on), tspan, x0, opts);

% Plot results
figure('Color','w','Position',[100 100 1400 600]);

subplot(1,2,1);
plot(t, x(:,1)*1000, 'b', 'LineWidth', 1.5); hold on;
plot(t, x(:,2)*1000, 'r', 'LineWidth', 1.5);
xline(t_sound_on, 'k--', 'LineWidth', 1.5, 'Label', 'Sound ON');
xlabel('Time [s]'); ylabel('Displacement [mm]');
legend('z (membrane)', 'y (frame)', 'Location','best');
title('Time History (Best Design: Rank 16)');
grid on;

subplot(1,2,2);
expansion = (x(:,2) - y_eq) * 1000;  % Frame expansion [mm]
plot(t, expansion, 'r', 'LineWidth', 2);
xline(t_sound_on, 'k--', 'LineWidth', 1.5, 'Label', 'Sound ON');
xlabel('Time [s]'); ylabel('Frame Expansion [mm]');
title('Expansion (y - y_{eq})');
grid on;

fprintf('\n========== SIMULATION RESULTS ==========\n');
fprintf('Design: x_push=%.3f, w_stiff=%.2f, w_mass=%.1f, L=%.1f\n', ...
    des_best.x_push_m*1000, des_best.w1_stiff_m*1000, ...
    des_best.w1_mass_m*1000, des_best.L_m*1000);
fprintf('Initial deflection z_eq = %.3f mm\n', z_eq*1000);
fprintf('Resonance frequency = %.1f Hz\n', p.freq);
fprintf('Max membrane amplitude = %.3f mm\n', max(x(:,1))*1000 - z_eq*1000);
fprintf('Max frame expansion = %.3f mm\n', max(expansion));
fprintf('=========================================\n');

function dx = eom_expansion(t, x, p, t_on)
    z = x(1); y = x(2); dz = x(3); dy = x(4);
    ct = y - (p.Gm * z^2);
    Fz = p.k1 * z - p.Kc * ct * (2 * p.Gm * z);
    Fy = p.k2 * (y - p.y_nat) + p.Kc * ct;
    if t >= t_on
        rp = min(1.0, (t - t_on) / 0.1);
        Fs = rp * p.P0 * (p.L * p.w1_mass) * cos(2*pi*p.freq*t);
    else
        Fs = 0;
    end
    ddz = (Fs - p.c1*dz - Fz) / p.m1;
    ddy = (0  - p.c2*dy - Fy) / p.m2;
    dx = [dz; dy; ddz; ddy];
end
