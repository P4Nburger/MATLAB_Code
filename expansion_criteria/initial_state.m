function [p, x0] = initial_state(p, x_push)
%INITIAL_STATE  押し込み量から釣り合い状態と初期条件を求める
%
%   [p, x0] = initial_state(p, x_push)
%     p      : MODEL_PARAMS で作ったパラメータ
%     x_push : フレームの押し込み量 [m]
%
%     p      : 釣り合い状態（x_push, z_eq, y_eq, y_natural）を追加したパラメータ
%     x0     : 初期状態 [z; y; dz/dt; dy/dt]（釣り合い点で静止）
%
%   初期たわみ z_eq は式3.5（x_push = C_eff*δ0^2/2、mm 単位の式）から求める。
%   フレームの自然長 y_natural は、z = z_eq, y = y_eq で力が釣り合うように決める。
%
%   See also MODEL_PARAMS, EOM_2DOF.

% 式3.5（mm 単位の式なので m と mm を換算する）
z_eq = sqrt(2 * (x_push * 1000) / p.C_eff) / 1000;

% z 方向の釣り合い: k1*z - 2*Γ*z*K_couple*(y - Γ*z^2) = 0
%   → y_eq - Γ*z_eq^2 = k1 / (2*Γ*K_couple)
coupling_eq = p.k1 / (2 * p.Gamma * p.K_couple);
y_eq = p.Gamma * z_eq^2 + coupling_eq;

% y 方向の釣り合い: k2*(y_eq - y_natural) + K_couple*(y_eq - Γ*z_eq^2) = 0
y_natural = y_eq + (p.K_couple * coupling_eq) / p.k2;

p.x_push    = x_push;
p.z_eq      = z_eq;
p.y_eq      = y_eq;
p.y_natural = y_natural;

x0 = [z_eq; y_eq; 0; 0];
end
