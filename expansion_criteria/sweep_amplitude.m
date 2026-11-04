function amp = sweep_amplitude(freqs, p, x0, tspan, ode_opts, use_parallel)
%SWEEP_AMPLITUDE  周波数ごとに膜のたわみの定常振幅を計算する
%
%   amp = sweep_amplitude(freqs, p, x0, tspan, ode_opts, use_parallel)
%     freqs        : 調べる周波数 [Hz]
%     p, x0        : INITIAL_STATE で作ったパラメータと初期状態
%     tspan        : 1周波数あたりの計算時間 [s]（音波は t = 0 から入れる）
%     ode_opts     : ode45 のオプション（[] なら ode45 の既定値）
%     use_parallel : true なら parfor で並列に計算する
%
%     amp          : 各周波数での振幅 [m]（計算時間の後半の (最大 - 最小)/2）
%
%   See also FIND_RESONANCE.

amp = zeros(size(freqs));
if use_parallel
    parfor i = 1:numel(freqs)
        amp(i) = steady_amplitude(freqs(i), p, x0, tspan, ode_opts);
    end
else
    for i = 1:numel(freqs)
        amp(i) = steady_amplitude(freqs(i), p, x0, tspan, ode_opts);
    end
end
end

function a = steady_amplitude(freq, p, x0, tspan, ode_opts)
p.freq = freq;
if isempty(ode_opts)
    [~, x] = ode45(@(t, x) eom_2dof(t, x, p, 0), tspan, x0);
else
    [~, x] = ode45(@(t, x) eom_2dof(t, x, p, 0), tspan, x0, ode_opts);
end
z_steady = x(round(end/2):end, 1);   % 後半を定常状態とみなす
a = (max(z_steady) - min(z_steady)) / 2;
end
