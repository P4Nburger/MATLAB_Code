# 音波駆動拡張ユニット構造の解析コード

明治大学 理工学部 機械情報工学科 材料システム研究室 卒業論文
**「音波駆動拡張ユニット構造の高性能化に向けた理論的検討」** で用いた MATLAB コードです。

> **English summary** — MATLAB code for the bachelor's thesis *"Theoretical study toward improving the performance of an acoustically driven expandable unit structure"* (Meiji University).
> It contains (1) a semi-empirical relation between the frame push-in and the initial deflection of the membrane, and (2) a two-degree-of-freedom model that predicts whether sound-wave excitation can trigger snap-through, including an automatic search for the critical push-in.

---

## フォルダ構成

```
.
├── initial_deflection_condition/  初期たわみ条件の算出（第3章）
├── expansion_criteria/            拡張可能条件の判定（第4章）
├── energy/                        弾性エネルギーの算出に用いた数値データ
├── biaxial_bending/               卒論提出後の追加検討（板の二軸曲げ・曲率評価）
├── tests/                         expansion_criteria/ のモデルのテスト
└── archive/                       開発途中の旧バージョン（参考用）
```

論文の結果に対応するのは `initial_deflection_condition/` と `expansion_criteria/` です。

---

## 論文との対応

| 内容 | スクリプト | 論文 |
|---|---|---|
| 押し込み量から初期たわみを求める | [`initial_deflection_condition/run_initial_deflection.m`](initial_deflection_condition/run_initial_deflection.m) | 第3章（式 3.5） |
| 1条件での周波数スイープと時間応答 | [`expansion_criteria/run_response_basic.m`](expansion_criteria/run_response_basic.m) | 第4章 |
| 上記の高速版（2段階スイープ・並列計算） | [`expansion_criteria/run_response_fast.m`](expansion_criteria/run_response_fast.m) | 第4章 |
| 拡張できる押し込み量の境界を自動探索 | [`expansion_criteria/run_boundary_search.m`](expansion_criteria/run_boundary_search.m) | 第4章 図4.2 |

論文を書いた時点のファイル名は、下の「[旧ファイル名との対応](#旧ファイル名との対応)」で調べられます。

---

## 実行方法

1. このリポジトリを clone するか、ZIP でダウンロードします。
2. MATLAB で実行したいスクリプトのあるフォルダに移動し、スクリプト名を入力して実行します。

```matlab
cd expansion_criteria
run_boundary_search
```

- 実行するのは `run_` で始まるスクリプトです。それ以外のファイルはスクリプトから呼ばれる関数です（[命名規則](#コードの構成と命名規則)）。
- 解析条件は各スクリプトの先頭の「1. 設定」で変更できます。構造・材料の値は [`expansion_criteria/model_params.m`](expansion_criteria/model_params.m) にまとめてあります。
- `run_response_fast.m` と `run_boundary_search.m` は、Parallel Computing Toolbox があれば並列計算します。なければ逐次計算になります（結果は同じで、時間が長くかかります）。

---

## 各フォルダの内容

### `initial_deflection_condition/` — 初期たわみ条件

FEA の結果で補正した理論式（式 3.5）を使い、フレームの押し込み量 $x_{\rm push}$ から膜の初期たわみ $\delta_0$ を求めます。

$$
x_{\rm push} = \frac{1}{2} C_{\rm eff}\, \delta_0^{2}
\quad\Longrightarrow\quad
\delta_0 = \sqrt{\frac{2\,x_{\rm push}}{C_{\rm eff}}}
\qquad [\mathrm{mm}]
$$

| ファイル | 説明 |
|---|---|
| `run_initial_deflection.m` | 押し込み量（1つ、またはベクトル）に対する初期たわみを計算する。モデルはひし形（$C_{\rm eff}=0.2795$）と短冊（$C_{\rm eff}=0.2548$）から選べる |

### `expansion_criteria/` — 拡張可能条件

膜の面外たわみ $z$ とフレームの変位 $y$ を自由度とする 2 自由度モデルで、音波を当てたときにスナップスルー（膜の反転）が起きるかを調べます。

**実行するスクリプト**

| ファイル | 説明 |
|---|---|
| `run_response_basic.m` | 基本版。周波数スイープで共振周波数を求め、その周波数で時間応答を計算する |
| `run_response_fast.m` | 高速版。1 Hz 刻みで広く調べてから、ピークの前後を細かく調べ直す（2段階スイープ） |
| `run_boundary_search.m` | 押し込み量を変えながら上の計算を繰り返し、スナップスルーが起きる押し込み量の境界を探索する（二分探索で範囲を絞り、刻み幅を 0.01 → 0.001 → 0.0001 mm と細かくする） |

**スクリプトから呼ばれる関数**

| ファイル | 説明 |
|---|---|
| `model_params.m` | 構造・材料パラメータの既定値と、そこから計算する質量・ばね定数・減衰係数 |
| `initial_state.m` | 押し込み量から釣り合い状態（初期たわみ・フレームの自然長）と初期条件を求める |
| `eom_2dof.m` | 運動方程式 |
| `sweep_amplitude.m` | 周波数ごとの定常振幅を計算する |
| `find_resonance.m` | 周波数スイープ（1段階または2段階）で共振周波数を求める |
| `simulate_response.m` | 指定した周波数の音波を当てたときの時間応答を計算する |
| `start_parallel_pool.m` | 並列計算の準備（Parallel Computing Toolbox がなければ逐次計算にする） |

<details>
<summary>モデルの概要（クリックで展開）</summary>

記号はコード中の変数名に合わせています。

$$
\begin{aligned}
m_1 \ddot z + c_1 \dot z + k_1 z - 2\Gamma K_{\rm c}\, z\,(y - \Gamma z^2) &= F(t) \\
m_2 \ddot y + c_2 \dot y + k_2 (y - y_{\rm nat}) + K_{\rm c}\,(y - \Gamma z^2) &= 0
\end{aligned}
$$

- 音波による力：$F(t) = r(t)\, P_0\, L\, w_{\rm mass} \cos(2\pi f t)$（$r(t)$ は音を入れてから 0.1 s で 0→1 になるランプ）
- $k_1 = 384 E_1 I_1 / L^3$、$I_1 = w_{\rm stiff} t_1^3 / 12$、$K_{\rm c} = E_1 w_{\rm stiff} t_1 / L$
- $m_1 = \rho_1 L w_{\rm mass} t_1$、$m_2 = \rho_2 S_2 L$、$c_i = 2\zeta\sqrt{m_i k_i}$
- $\Gamma = C_{\rm eff}/2$（式 3.5 の $x_{\rm push}$–$\delta_0$ 関係に対応）
- 初期状態は $z = \delta_0$ で静止した釣り合い点。釣り合い条件からフレームの自然長 $y_{\rm nat}$ を決める
- スナップスルーの判定：時間応答で $\min z(t) < -0.8\, \delta_0$ となれば「あり」

</details>

主なパラメータ（`model_params.m` の既定値）：

| 記号 | 変数名 | 値 |
|---|---|---|
| 膜の長さ $L$ | `L` | 18 mm |
| 剛性の計算に使う幅（最も細い部分） | `w1_stiff` | 5.0 mm |
| 質量・受音面積の計算に使う幅 | `w1_mass` | 22.5 mm |
| 膜厚 $t_1$ | `t1` | 0.1 mm |
| ヤング率 $E_1$ | `E1` | 3.45 GPa |
| 密度 $\rho_1, \rho_2$ | `rho1`, `rho2` | 1250 kg/m³ |
| フレームの断面積 $S_2$ | `S2` | 152.522 mm² |
| フレームのばね定数 $k_2$ | `k2` | 222 N/m |
| 減衰比 $\zeta$ | `zeta` | 0.05 |
| 有効形状定数 $C_{\rm eff}$ | `C_eff` | 0.2794 |
| 音圧振幅 $P_0$ | `P0` | 2 Pa |

値を変えるときは、`model_params('P0', 5)` のように名前と値の組で上書きします。

### `energy/` — 弾性エネルギー算出用データ

弾性エネルギーの算出に用いた静解析の結果（2 列の数値データ）です。

| ファイル | 説明 |
|---|---|
| `static_hole.csv` | 穴ありモデル |
| `static_normal.csv` | 通常モデル |

### `biaxial_bending/` — 卒論提出後の追加検討

卒論の範囲外です。膜を板として扱った二軸曲げと、根本の曲率の評価を検討しています。

| ファイル | 説明 |
|---|---|
| `plate_deflection_levy.m` | x=0 固定・x=lx 自由の板に等分布荷重をかけたときのたわみ分布（y 方向に sin 展開した解の第 1 項） |
| `plate_deflection_levy_uniform_y.m` | 上と同じ x 方向の解を y 方向に一様に広げた版 |
| `derive_deflection_coefficients.m` | たわみ式の係数 A〜D の記号計算による導出 |
| `plate_deflection_ritz.m` | Ritz 法による板のたわみ・曲げモーメント分布 |
| `curvature_from_3points.m` | 3 点を通る円から曲率半径・曲率を求める |
| `root_curvature_cubic_fit.m` | 根本から 6 点を 3 次式でフィットして根本の曲率を求める（片持ちはりと部分固定の比較） |

### `tests/` — テスト

`expansion_criteria/` のモデルについて、釣り合い点で加速度が 0 になること、初期たわみが式 3.5 を満たすことなどを確かめます。リポジトリのルートで次を実行します。

```matlab
results = runtests('tests')
```

### `archive/` — 旧バージョン

開発の途中で作ったコードです。論文の結果には使っておらず、保守もしていません。詳しくは [`archive/README.md`](archive/README.md) を見てください。

---

## コードの構成と命名規則

- **ファイル名は英語の snake_case** で、中身がわかる名前にしています（例：`run_boundary_search.m`、`initial_state.m`）。
- **`run_` で始まるファイルが実行用のスクリプト**です。それ以外は、スクリプトから呼ばれる関数です。
- **モデルの定義は1か所にまとめています。** パラメータは `model_params.m`、運動方程式は `eom_2dof.m` にあります。どのスクリプトも同じものを使うので、スクリプトごとに値がずれることはありません。
- **単位**：`expansion_criteria/` は SI 単位（m, kg, s）、`initial_deflection_condition/` と `biaxial_bending/` の曲率計算は mm です。
- 各ファイルの先頭に説明を書いています。MATLAB で `help run_boundary_search` のように実行すると表示されます。

---

## 旧ファイル名との対応

論文を書いた時点のファイル名と、現在のファイル名の対応です。
**名前の整理（2026/11）で計算内容は変えていません。** 旧コードと新コードで同じ計算（計算時間を短くした条件、GNU Octave で実行）を行い、周波数スイープの結果・時間応答・境界探索の手順がビット単位で一致することを確認しています。

| 旧ファイル名 | 現在のファイル |
|---|---|
| `initial_deflection_condition/current/kotei_LikoruS_danmen.m` | `initial_deflection_condition/run_initial_deflection.m` |
| `initial_deflection_condition/current/kotei_LikoruS_danmen2.m` | `initial_deflection_condition/run_initial_deflection.m`（`C_eff = 0.2781`、`x_push = 0.01:0.01:0.1` に設定） |
| `expansion_criteria/current/gouseihenkou.m` | `expansion_criteria/run_response_basic.m` |
| `expansion_criteria/current/機能追加版/scan.m` | `expansion_criteria/run_response_fast.m` |
| `expansion_criteria/current/機能追加版/run_simulation_optimized.m` | `expansion_criteria/run_boundary_search.m` |
| `expansion_criteria/current/機能追加版/equations_2DOF_SelfSustain.m` | `expansion_criteria/eom_2dof.m` |
| `expansion_criteria/current/機能追加版/run_simulation.m` | `archive/expansion_criteria/v8/run_simulation.m` |
| `expansion_criteria/current/onatuhenkou.m` | `archive/expansion_criteria/after_thesis/onatuhenkou.m`（`run_response_fast.m` で `P0 = 5`、`x_push = 0.0735e-3`、`refine_step = 0.1` にすると同じ計算） |
| `2軸曲げ/nijikumage.m` | `biaxial_bending/plate_deflection_levy.m` |
| `2軸曲げ/nijikumageynasi.m` | `biaxial_bending/plate_deflection_levy_uniform_y.m` |
| `2軸曲げ/seiteki.m` | `biaxial_bending/derive_deflection_coefficients.m` |
| `2軸曲げ/kinji.m` | `biaxial_bending/plate_deflection_ritz.m` |
| `2軸曲げ/kyokuritukeisan.m` | `biaxial_bending/curvature_from_3points.m` |
| `2軸曲げ/current/kyokuritu4ten.m` | `biaxial_bending/root_curvature_cubic_fit.m` |
| `エネルギー/` | `energy/` |

パラメータの変数名も一部変えています：`delta`（減衰比）→ `zeta`、`k2_base` → `k2`。

---

## 動作環境

- **MATLAB** R2025b（筆者の環境）。R2019b 以降なら動く想定です（スクリプト内のローカル関数、`yline`、`tiledlayout` を使っているため）
- **Parallel Computing Toolbox**（任意）：あれば `run_response_fast.m`、`run_boundary_search.m` が並列計算になります
- **Symbolic Math Toolbox**：`biaxial_bending/derive_deflection_coefficients.m` で必要
- `archive/` の一部は Optimization Toolbox・Curve Fitting Toolbox を使います

---

## 更新履歴

- **2026/1/1**：初版
- **2026/1/29**：コード追加、ファイル整理、README 修正
- **2026/1/31**：ファイル名の整理、README 修正
- **2026/3/4**：`energy/` を追加
- **2026/3〜4**：卒論提出後の追加検討（音圧 5 Pa の検討、二軸曲げ）を追加
- **2026/11/4**：リポジトリを整理
  - フォルダ構成の整理（`current/` を廃止、旧版を `archive/` に集約、フォルダ名を英語に統一、MATLAB プロジェクトの自動生成ファイルを削除）
  - ファイル名を英語の snake_case に統一し、`run_` で始まるスクリプトと関数を区別
  - `expansion_criteria/` で重複していたパラメータと運動方程式を共通の関数にまとめた（計算結果は変わらないことを確認済み）
  - Parallel Computing Toolbox がなくても動くようにした
  - `run_initial_deflection.m` で短冊モデルを選ぶとエラーになっていた問題を修正
  - テスト（`tests/`）を追加、README を全面改訂

---

## Author

P4Nburger
Meiji University — Material System Lab
