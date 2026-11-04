function [t, x] = simulate_response(p, x0, freq, tspan, t_sound_on, ode_opts)
%SIMULATE_RESPONSE  指定した周波数の音波を当てたときの時間応答を計算する
%
%   [t, x] = simulate_response(p, x0, freq, tspan, t_sound_on, ode_opts)
%     p, x0      : INITIAL_STATE で作ったパラメータと初期状態
%     freq       : 音波の周波数 [Hz]
%     tspan      : 計算時間 [s]
%     t_sound_on : 音波を入れ始める時刻 [s]
%     ode_opts   : ode45 のオプション
%
%     t, x       : 時刻 [s] と状態 [z, y, dz/dt, dy/dt]（各行が1時刻）
%
%   See also EOM_2DOF.

p.freq = freq;
[t, x] = ode45(@(t, x) eom_2dof(t, x, p, t_sound_on), tspan, x0, ode_opts);
end
