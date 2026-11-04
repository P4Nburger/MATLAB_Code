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
└── archive/                       開発途中の旧バージョン（参考用）
```

論文の結果に対応するのは `initial_deflection_condition/` と `expansion_criteria/` です。
`archive/` は開発の経緯を残すためのもので、保守していません。

---

## 論文との対応

| 内容 | スクリプト | 論文 |
|---|---|---|
| 押し込み量から初期たわみを求める | [`initial_deflection_condition/kotei_LikoruS_danmen.m`](initial_deflection_condition/kotei_LikoruS_danmen.m) | 第3章（式 3.5） |
| 1条件での周波数スイープと時間応答 | [`expansion_criteria/gouseihenkou.m`](expansion_criteria/gouseihenkou.m) | 第4章 |
| 上記の高速版（並列計算・2段階スイープ） | [`expansion_criteria/scan.m`](expansion_criteria/scan.m) | 第4章 |
| 拡張できる押し込み量の境界を自動探索 | [`expansion_criteria/run_simulation_optimized.m`](expansion_criteria/run_simulation_optimized.m) | 第4章 図4.2 |

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
| `kotei_LikoruS_danmen.m` | ひし形モデル（$C_{\rm eff}=0.2795$）で、指定した押し込み量に対する初期たわみを計算 |
| `kotei_LikoruS_danmen2.m` | $C_{\rm eff}=0.2781$ を使い、押し込み量 0.01〜0.1 mm の範囲で計算する版 |

### `expansion_criteria/` — 拡張可能条件

膜の面外たわみ $z$ とフレームの変位 $y$ を自由度とする 2 自由度モデルで、音波を当てたときにスナップスルー（膜の反転）が起きるかを調べます。

| ファイル | 説明 |
|---|---|
| `gouseihenkou.m` | 基本コード。周波数スイープで共振周波数を求め、その周波数で時間応答を計算する |
| `scan.m` | `gouseihenkou.m` の高速版。粗いスイープ（1 Hz 刻み）→ ピーク周辺の細かいスイープの 2 段階で共振周波数を求める |
| `run_simulation_optimized.m` | 押し込み量を変えながら `scan.m` と同様の計算を繰り返し、スナップスルーが起きる押し込み量の境界を探索する（二分探索で範囲を絞り、刻み幅を段階的に細かくする） |
| `equations_2DOF_SelfSustain.m` | 運動方程式（上の各スクリプトにも同じ関数がローカル関数として入っている） |
| `onatuhenkou.m` | 卒論提出後に追加。音圧を 5 Pa にした場合の検討 |

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

主なパラメータ（コード内の値）：

| 記号 | 変数名 | 値 |
|---|---|---|
| 膜の長さ $L$ | `L` | 18 mm |
| 剛性計算用の幅（最も細い部分） | `w1_stiff` | 5.0 mm |
| 質量・受音面積用の幅 | `w1_mass` | 22.5 mm |
| 膜厚 $t_1$ | `t1` | 0.1 mm |
| ヤング率 $E_1$ | `E1` | 3.45 GPa |
| 密度 $\rho_1, \rho_2$ | `rho1`, `rho2` | 1250 kg/m³ |
| フレーム断面積 $S_2$ | `S2` | 152.522 mm² |
| フレームのばね定数 $k_2$ | `k2_base` | 222 N/m |
| 減衰比 $\zeta$ | `delta` | 0.05 |
| 音圧振幅 $P_0$ | `P0_drive` | 2 Pa |
| 有効形状定数 $C_{\rm eff}$ | `C_eff` | 0.2794 |

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
| `nijikumage.m` | x=0 固定・x=lx 自由の板に等分布荷重をかけたときのたわみ分布（y 方向に sin 展開した解の第 1 項） |
| `nijikumageynasi.m` | 上と同じ x 方向の解を y 方向に一様に広げた版 |
| `seiteki.m` | たわみ式の係数 A〜D の記号計算による導出 |
| `kinji.m` | Ritz 法による板のたわみ・曲げモーメント分布 |
| `kyokuritukeisan.m` | 3 点を通る円から曲率半径・曲率を求める |
| `kyokuritu4ten.m` | 根本から 6 点を 3 次式でフィットして根本の曲率を求める（片持ちはりと部分固定の比較） |

### `archive/` — 旧バージョン

開発の途中で作ったコードです。詳しくは [`archive/README.md`](archive/README.md) を見てください。

---

## 実行方法

1. このリポジトリを clone するか、ZIP でダウンロードします。
2. MATLAB で実行したいスクリプトのあるフォルダに移動し、スクリプト名を入力して実行します。

```matlab
cd expansion_criteria
run_simulation_optimized
```

- 各スクリプトは単体で動きます（必要な関数はファイル内にローカル関数として入っています）。パスの設定は不要です。
- パラメータは各スクリプトの先頭の「パラメータ設定」で変更できます。
- `scan.m`、`run_simulation_optimized.m`、`onatuhenkou.m` は並列計算（`parpool` / `parfor`）を使うため、Parallel Computing Toolbox が必要です。特に `run_simulation_optimized.m` は計算量が多いので、並列計算できる環境で実行してください。

---

## 動作環境

- **MATLAB** R2025b（筆者の環境）
- **Parallel Computing Toolbox**：`scan.m`、`run_simulation_optimized.m`、`onatuhenkou.m` で必要
- **Symbolic Math Toolbox**：`biaxial_bending/seiteki.m` で必要
- `archive/` の一部は Optimization Toolbox・Curve Fitting Toolbox を使います

---

## 更新履歴

- **2026/1/1**：初版
- **2026/1/29**：コード追加、ファイル整理、README 修正
- **2026/1/31**：ファイル名の整理、README 修正
- **2026/3/4**：`energy/` を追加
- **2026/3〜4**：卒論提出後の追加検討（`onatuhenkou.m`、`biaxial_bending/`）を追加
- **2026/11/4**：リポジトリ構成を整理（`current/` を廃止して最新版を各フォルダ直下に、旧版を `archive/` に集約、フォルダ名を英語に統一、MATLAB プロジェクトの自動生成ファイルを削除）、README を全面改訂

---

## Author

P4Nburger
Meiji University — Material System Lab
