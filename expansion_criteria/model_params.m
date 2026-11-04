function p = model_params(varargin)
%MODEL_PARAMS  2自由度モデルのパラメータ（既定値）と派生定数を返す
%
%   p = model_params()                 既定値を使う
%   p = model_params('P0', 5, ...)     名前と値の組で既定値を上書きする
%
%   単位は SI（m, kg, s, N, Pa）。
%   既定値は論文 第4章の解析条件。派生定数（質量・ばね定数・減衰係数など）は
%   上書きした後の値から計算する。
%
%   See also INITIAL_STATE, EOM_2DOF.

% --- 既定値 ---
p = struct( ...
    'L',        18e-3, ...       % [m]      膜の長さ
    'w1_stiff', 5.0e-3, ...      % [m]      剛性の計算に使う幅（最も細い部分）
    'w1_mass',  22.5e-3, ...     % [m]      質量・受音面積の計算に使う幅（平均幅）
    't1',       0.1e-3, ...      % [m]      膜厚
    'E1',       3.45e9, ...      % [Pa]     膜のヤング率
    'rho1',     1250, ...        % [kg/m^3] 膜の密度
    'S2',       152.522e-6, ...  % [m^2]    フレームの断面積
    'rho2',     1250, ...        % [kg/m^3] フレームの密度
    'k2',       222.0, ...       % [N/m]    フレームのばね定数
    'zeta',     0.05, ...        % [-]      減衰比
    'C_eff',    0.2794, ...      % [-]      有効形状定数（式3.5）
    'P0',       2, ...           % [Pa]     音圧振幅
    'freq',     0);              % [Hz]     駆動周波数（解析の中で設定する）

% --- 上書き ---
if mod(numel(varargin), 2) ~= 0
    error('model_params:args', '引数は名前と値の組で指定してください。');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~isfield(p, name)
        error('model_params:unknown', '不明なパラメータ名です: %s', name);
    end
    p.(name) = varargin{k+1};
end

% --- 派生定数 ---
p.m1       = p.rho1 * p.L * p.w1_mass * p.t1;    % [kg]  膜の質量
p.m2       = p.rho2 * p.S2 * p.L;                % [kg]  フレームの質量
p.I1       = p.w1_stiff * p.t1^3 / 12;           % [m^4] 膜の断面二次モーメント
p.k1       = 384 * p.E1 * p.I1 / (p.L^3);        % [N/m] 膜の曲げのばね定数
p.K_couple = p.E1 * (p.w1_stiff * p.t1) / p.L;   % [N/m] 膜の伸びのばね定数（連成項）
p.c1       = 2 * p.zeta * sqrt(p.m1 * p.k1);     % [N s/m] 膜の減衰係数
p.c2       = 2 * p.zeta * sqrt(p.m2 * p.k2);     % [N s/m] フレームの減衰係数
p.Gamma    = (p.C_eff / 2) * 1000;               % [1/m] 幾何学係数（y ≈ Γ z^2、式3.5 を m 単位に換算）
end
