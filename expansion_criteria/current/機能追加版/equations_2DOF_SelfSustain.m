function dx = equations_2DOF_SelfSustain(t, x, p, t_sound_on)
z  = x(1);
y  = x(2);
dz = x(3);
dy = x(4);

k1_curr = p.k1;

coupling_term = y - (p.Gamma * z^2);

F_restore_z = k1_curr * z ...
    - p.K_couple * coupling_term * (2 * p.Gamma * z);

F_restore_y = p.k2_base * (y - p.y_natural) ...
    + p.K_couple * coupling_term;

if t >= t_sound_on
    ramp = min(1.0, (t - t_sound_on)/0.1);
    F_sound = ramp * p.P0 * (p.L * p.w1_mass) * cos(2 * pi * p.freq * t);
else
    F_sound = 0;
end

ddz = (F_sound - p.c1 * dz - F_restore_z) / p.m1;
ddy = (- p.c2 * dy - F_restore_y) / p.m2;

dx = [dz; dy; ddz; ddy];
end