% ==========================================
% 膜のたわみ式 係数決定のシンボリック導出
% ==========================================
syms E nu h lx ly qn alpha D a x A B C D_coeff

% 1. 行列 M の定義 (資料9ページ) [cite: 207]
M = [1, 0, 1, 0;
     0, 1, 0, 1;
     (1+nu)*cosh(a), (1+nu)*sinh(a), (nu-1)*cos(a), (nu-1)*sin(a);
     (3-nu)*sinh(a), (3-nu)*cosh(a), (nu-1)*sin(a), (1-nu)*cos(a)];

% 2. 右辺ベクトル b の定義 [cite: 205-206]
% ※ qn / (D * alpha^4) を便宜上 K と置くと見やすくなります
syms K
b = [-K; 0; -nu*K; 0];

% 3. 連立方程式を解く (係数 A, B, C, D の算出)
coeffs = M \ b;

% 4. 結果の表示
fprintf('--- 係数 A, B, C, D の導出結果 ---\n');
disp('A = '); disp(simplify(coeffs(1)));
disp('B = '); disp(simplify(coeffs(2)));
disp('C = '); disp(simplify(coeffs(3)));
disp('D = '); disp(simplify(coeffs(4)));

% 5. たわみ関数 Wn(x) の一般形を表示 [cite: 179]
Wn_x = coeffs(1)*cosh(alpha*x) + coeffs(2)*sinh(alpha*x) + ...
       coeffs(3)*cos(alpha*x) + coeffs(4)*sin(alpha*x) + K;
fprintf('\n--- たわみ関数 Wn(x) の一般式 ---\n');
disp(simplify(Wn_x));