% -------------------------
% 1. パラメータ設定
% (変更なし)
% -------------------------
clear; clc;

z0_initial   = 2e-3;
dz0_initial  = 0;
tspan = [0 20];
freq_range = 10:1:200;

alpha = 0.8;
gamma = 0.3;
beta  = 500;

params = struct( ...
        'L', 30e-3, 'P0', 2, 'E1', 197e9, 'w1', 10e-3, 't1', 0.02e-3, ...
        'roh1', 7.93e3, 'delta1', 0.02, 'w2', 7.5e-3, 'k2', 177.035, ...
        'delta2', 0.0278, 'L_side', 60e-3, 'L_up', 15e-3, 't2', 1e-3, ...
        'S2', 126.87e-6, 'roh2', 1.24e3, ...
        'alpha', alpha, 'gamma', gamma, 'beta', beta ...
    );

x0 = [z0_initial; dz0_initial; 0];

% -------------------------
% 2. シミュレーション実行とプロット
% (変更なし)
% -------------------------
params = generate_params(params);

z_init = z0_initial;
denom = params.L^2 - 2*z_init^2;
restoring_force_at_z0 = - (params.k1*z_init + (4*params.k2/denom)*z_init^3);
params.F_static_initial = -restoring_force_at_z0;

amp = freq_sweep(params, freq_range, tspan, x0);
[~, idx] = max(amp);
f_res = freq_range(idx);
fprintf('共振周波数: %.2f Hz\n', f_res);

omega_res = 2*pi*f_res;
[t, x] = ode45(@(t,x) nonlinearForcedODE(t,x,params,omega_res), tspan, x0);
z_sim = x(:,1);
A_sim = x(:,3);

% (時間応答グラフのプロット部分は変更なし)
figure;
plot(t, z_sim*1000, 'b-');
hold on;
center_of_vibration = z0_initial * exp(-params.beta * A_sim);
plot(t, center_of_vibration*1000, 'r--', 'LineWidth', 2.5);
yline(0, 'k--');
hold off;
set(gca,'Fontsize',16,'FontWeight','bold');
xlabel('時間 [s]'); ylabel('変位 z(t) [mm]');
title('【最終版】振動の大きさに応じて振動中心が移動する様子');
legend('z(t)の変位', '振動中心の移動', '最終的な振動中心');
grid on;


% -------------------------------------------------------------------
% ★★★ ここからが修正したアニメーションのコード ★★★
% -------------------------------------------------------------------

% y方向変位の計算
y_sim = -2*z_sim.^2 / params.L;

% アニメーションの準備
figure; % 新しいウィンドウ
% ★ z軸とy軸を入れ替えてプロット
h = plot(z_sim(1)*1000, y_sim(1)*1000, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
grid on;
hold on;

% ★ ラベルとタイトルを修正
xlabel('z変位 [mm]');
ylabel('y変位 [mm]');
title('膜中央の z-y 動き (アニメーション)');

% 表示範囲を固定する
max_z_display = max(abs(z_sim)) * 1.2 * 1000;
max_y_display = max(abs(y_sim)) * 1.2 * 1000;
% ★ xlim と ylim を入れ替え
xlim([-max_z_display, max_z_display]);
ylim([-max_y_display, 0.1*max_y_display]); % yは常に負なので、表示範囲を調整

axis equal; % 軸のスケールを同じに
plot(0, 0, 'k+', 'MarkerSize', 10); % 原点(0,0)に十字マーク
line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); % z=0 補助線
line(xlim, [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); % y=0 補助線
hold off;

% ★★★ アニメーションループをより安定な方式に変更 ★★★
fprintf('アニメーションを開始します...\n');
% 全ての点を描画すると重いので、間引いて表示する
step = 50; % 50点ごとに1フレームを描画
for i = 1:step:length(t)
    % ★ set の XData と YData を入れ替え
    set(h, 'XData', z_sim(i)*1000, 'YData', y_sim(i)*1000);
    
    % フレームを描画 (より安定な 'drawnow' を使用)
    drawnow;
    
    % 一定時間待つ (PCの性能に依存しないように)
    pause(0.01); 
end
fprintf('アニメーションが終了しました。\n');
% ★★★ ここまでが修正したアニメーションのコード ★★★

% -----------------------------------------------------
% --- ここから下は関数の定義 (変更なし) ---
% -----------------------------------------------------

function p = generate_params(p)
    p.I1 = p.w1 * p.t1^3 / 12; p.k1 = 384 * p.E1 * p.I1 / (p.L^3);
    p.m1 = p.roh1 * p.L * p.w1 * p.t1; p.c1 = p.delta1 / pi * sqrt(p.m1 * p.k1);
    p.S2 = 126.87e-6 * ((p.L_side + p.L_up)/75e-3) * (p.L/30e-3);
    p.m2 = p.S2 * p.w2 * p.roh2; p.c2 = p.delta2 / pi * sqrt(p.m2 * p.k2);
end

function dx = nonlinearForcedODE(t,x,p,omega)
    z = x(1);
    dz = x(2);
    A = x(3);
    
    F_static = p.F_static_initial * exp(-p.beta * A);
    F_dynamic = p.P0 * (p.L - 2* p.t2) * p.w1 * cos(omega*t);
    F_total = F_dynamic + F_static;
    
    denom = p.L^2 - 2*z^2;
    ddz = (1 / (p.m1 + 8*p.m2*z^2/denom)) * ...
         ( - (c1 + 8*p.c2*z^2/denom)*dz - p.k1*z ...
           - (8*p.m2*z/denom)*dz^2 - (4*p.k2/denom)*z^3 + F_total );
           
    dA = p.alpha * abs(z) - p.gamma * A;
    
    dx = [dz; ddz; dA];
end

function amp = freq_sweep(p, freq_range, tspan, x0)
    amp = zeros(size(freq_range));
    for i = 1:length(freq_range)
        omega = 2*pi*freq_range(i);
        [~,x] = ode45(@(t,x) nonlinearForcedODE(t,x,p,omega), tspan, x0);
        z = x(:,1);
        amp(i) = max(abs(z(round(end/2):end)));
    end
end