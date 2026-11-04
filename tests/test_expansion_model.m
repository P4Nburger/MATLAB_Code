function tests = test_expansion_model
%TEST_EXPANSION_MODEL  2自由度モデル（expansion_criteria/）の基本的な性質を確かめるテスト
%
%   実行方法（リポジトリのルートで）:
%       results = runtests('tests')

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
    fullfile(here, '..', 'expansion_criteria')));
end

function testDerivedConstants(testCase)
% 既定値から派生定数が計算されている
p = model_params();
verifyEqual(testCase, p.k1, 384 * p.E1 * (p.w1_stiff * p.t1^3 / 12) / p.L^3, 'RelTol', 1e-12);
verifyEqual(testCase, p.Gamma, p.C_eff / 2 * 1000, 'RelTol', 1e-12);
end

function testOverride(testCase)
% 名前と値の組で上書きでき、派生定数も上書き後の値で計算される
p = model_params('P0', 5, 't1', 0.2e-3);
verifyEqual(testCase, p.P0, 5);
verifyEqual(testCase, p.m1, p.rho1 * p.L * p.w1_mass * 0.2e-3, 'RelTol', 1e-12);
verifyError(testCase, @() model_params('no_such_param', 1), 'model_params:unknown');
end

function testInitialDeflectionFollowsEq35(testCase)
% 初期たわみが式3.5（x_push = C_eff*δ0^2/2、mm 単位）を満たす
x_push = 0.07e-3;
[p, x0] = initial_state(model_params(), x_push);
delta0_mm = p.z_eq * 1000;
verifyEqual(testCase, 0.5 * p.C_eff * delta0_mm^2, x_push * 1000, 'RelTol', 1e-12);
verifyEqual(testCase, x0, [p.z_eq; p.y_eq; 0; 0]);
end

function testEquilibrium(testCase)
% 音波がなければ、初期状態（釣り合い点）で加速度が 0
[p, x0] = initial_state(model_params(), 0.07e-3);
dx = eom_2dof(0, x0, p, Inf);
scale = p.k1 * p.z_eq / p.m1;   % 加速度の大きさの目安 [m/s^2]
verifyLessThan(testCase, max(abs(dx)), 1e-9 * scale);
end

function testStaysAtRestWithoutSound(testCase)
% 音波がなければ、釣り合い点から動かない
[p, x0] = initial_state(model_params(), 0.07e-3);
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-12);
[~, x] = simulate_response(p, x0, 30, [0 0.2], Inf, opts);
verifyLessThan(testCase, max(abs(x(:, 1) - p.z_eq)), 1e-6 * p.z_eq);
end
