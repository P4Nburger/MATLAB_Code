function dx = eom_2dof(t, x, p, t_sound_on)
%EOM_2DOF  2自由度モデルの運動方程式（ode45 に渡す関数）
%
%   dx = eom_2dof(t, x, p, t_sound_on)
%     x          : 状態 [z; y; dz/dt; dy/dt]
%                  z = 膜のたわみ [m]、y = フレームの変位 [m]
%     p          : MODEL_PARAMS と INITIAL_STATE で作ったパラメータ
%     t_sound_on : 音波を入れ始める時刻 [s]（0.1 s かけて振幅を立ち上げる）
%
%   運動方程式
%     m1*z'' + c1*z' + k1*z - 2*Γ*K_couple*z*(y - Γ*z^2) = F(t)
%     m2*y'' + c2*y' + k2*(y - y_natural) + K_couple*(y - Γ*z^2) = 0
%     F(t) = r(t) * P0 * L * w1_mass * cos(2*pi*f*t)   （r(t): 立ち上げのランプ）
%
%   See also MODEL_PARAMS, INITIAL_STATE, SIMULATE_RESPONSE.

T_RAMP = 0.1;   % [s] 音波の振幅を 0 から 1 まで上げる時間

z  = x(1);
y  = x(2);
dz = x(3);
dy = x(4);

% 復元力
coupling_term = y - (p.Gamma * z^2);
F_restore_z = p.k1 * z - p.K_couple * coupling_term * (2 * p.Gamma * z);
F_restore_y = p.k2 * (y - p.y_natural) + p.K_couple * coupling_term;

% 音波による外力（音を受ける面積は平均幅 w1_mass で計算する）
if t >= t_sound_on
    ramp    = min(1.0, (t - t_sound_on) / T_RAMP);
    F_sound = ramp * p.P0 * (p.L * p.w1_mass) * cos(2 * pi * p.freq * t);
else
    F_sound = 0;
end

ddz = (F_sound - p.c1 * dz - F_restore_z) / p.m1;
ddy = (0 - p.c2 * dy - F_restore_y) / p.m2;

dx = [dz; dy; ddz; ddy];
end
