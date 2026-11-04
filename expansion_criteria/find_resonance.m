function res = find_resonance(p, x0, sweep)
%FIND_RESONANCE  周波数スイープで共振周波数を求める
%
%   res = find_resonance(p, x0, sweep)
%
%   1段目で sweep.freqs を調べて振幅が最大の周波数を探し、
%   2段目でその前後 ±refine_range を refine_step 刻みで調べ直す。
%
%   sweep のフィールド
%     freqs        : 1段目に調べる周波数 [Hz]
%     refine_range : 2段目でピークの前後に調べる幅 [Hz]（[] なら2段目を省略）
%     refine_step  : 2段目の刻み幅 [Hz]
%     tspan        : 1周波数あたりの計算時間 [s]
%     ode_opts     : ode45 のオプション（[] なら ode45 の既定値）
%     use_parallel : true なら並列計算
%
%   res のフィールド
%     f_res, amp_res     : 共振周波数 [Hz] とその振幅 [m]
%     freqs1, amp1       : 1段目の周波数と振幅
%     f_peak1, amp_peak1 : 1段目のピーク
%     freqs2, amp2       : 2段目の周波数と振幅（省略したときは空）
%
%   See also SWEEP_AMPLITUDE.

% 1段目
res.freqs1 = sweep.freqs;
res.amp1   = sweep_amplitude(res.freqs1, p, x0, sweep.tspan, sweep.ode_opts, sweep.use_parallel);
[res.amp_peak1, idx] = max(res.amp1);
res.f_peak1 = res.freqs1(idx);

res.f_res   = res.f_peak1;
res.amp_res = res.amp_peak1;
res.freqs2  = [];
res.amp2    = [];
if isempty(sweep.refine_range)
    return
end

% 2段目（1段目の範囲からははみ出さない）
f_start    = max(min(sweep.freqs), res.f_peak1 - sweep.refine_range);
f_end      = min(max(sweep.freqs), res.f_peak1 + sweep.refine_range);
res.freqs2 = f_start:sweep.refine_step:f_end;
res.amp2   = sweep_amplitude(res.freqs2, p, x0, sweep.tspan, sweep.ode_opts, sweep.use_parallel);
[res.amp_res, idx] = max(res.amp2);
res.f_res = res.freqs2(idx);
end
