# Acoustic Metamaterial Analysis Code
明治大学 理工学部 機械情報工学科 材料システム研究室 卒業論文
「音波駆動拡張ユニット構造の高性能化に向けた理論的検討」解析コードリポジトリ

## 概要
本リポジトリは、音波駆動型拡張メタマテリアル構造の設計・解析を行うための MATLAB コード群です。
主に以下の2つの解析を行います。

1. **初期たわみ条件の算出**: 幾何学的非線形性を考慮した理論式に基づく計算（第3章対応）
2. **拡張可能条件の判定**: 2自由度連成振動モデルを用いたスナップスルー挙動の数値解析（第4章対応）

---

## フォルダ構成

### 📂 `initial_deflection_condition/`
**初期たわみ条件を求めるコード群**
押し込み量と初期たわみ量の関係を解析します。

- **`current/`**: 最新版コード（実行推奨）
- **`archive/`**: 過去バージョン・実験的コード

#### 主要スクリプト（current/）
| ファイル名 | 説明 |
|-----------|------|
| `kotei_LikoruS_danmen.m` | **初期たわみ算出コード**<br>有限要素解析で補正した理論式（式 3.5）を用い、任意の押し込み量に対する初期たわみを計算します。 |

---

### 📂 `expansion_criteria/`
**拡張可能条件を判定するコード群**
音波照射下でのスナップスルー発生有無を判定し、拡張可能な初期たわみ領域を探索します。

- **`current/`**: 最新版コード
  - **`機能追加版/`**: 並列計算対応・最適化版（**推奨**）
- **`archive/`**: 過去バージョン・実験的コード

#### 基本解析スクリプト（current/）
| ファイル名 | 説明 |
|-----------|------|
| `gouseihenkou.m` | 拡張判定の基本コード（単一パラメータ確認用） |

#### 機能追加版解析スクリプト（current/機能追加版/）
| ファイル名 | 説明 |
|-----------|------|
| `scan.m` | **高速判定ソルバー**<br>`gouseihenkou.m` の計算アルゴリズムを最適化したもの。高速かつ高精度に判定を行います。 |
| `run_simulation_optimized.m` | **拡張可能範囲の自動探索**<br>並列計算を用いてパラメータスイープを行い、スナップスルー可能な押し込み量の境界値を探索します。（論文 図4.3に対応） |
| `run_simulation.m` | `run_simulation_optimized.m` のプロトタイプ版 |

---

## 使い方（推奨フロー）

### 1. 初期たわみ条件の計算
```matlab
cd initial_deflection_condition/current/
kotei_LikoruS_danmen  
```

### 2. 2. 拡張可能条件の探索（詳細解析）
```matlab
cd expansion_criteria/current/機能追加版/
run_simulation_optimized  % 並列計算による高速探索
```

### 3. 拡張可能な初期たわみ量の詳細な探索
```matlab
cd expansion_criteria/current/
scan  % 任意の初期たわみでの解析実行
```

---

## 環境

- **MATLAB**: R2025bを推奨(筆者の環境)
- **必須 Toolbox**: 
  - Optimization Toolbox
- **推奨 Toolbox**:
  - Parallel Computing Toolbox
> Note: run_simulation_optimized.m などの探索コードは計算負荷が高いため、Parallel Computing Toolbox による並列計算環境での実行を強く推奨します。
---

## 注意事項

- `current/` フォルダ内が最新版です。`archive/` は開発履歴の保存用です。
- 論文中の解析結果再現には、並列計算環境の使用を前提としたパラメータ設定となっている箇所があります。



---

## 更新履歴

- **2026/1/1**: 初期バージョン作成
- **2026/1/29**: コード追加　ファイル整理、README 修正
- **2026/1/31**: ファイル名前整理、README 修正
---

## Author

P4Nburger  
Meiji University - Material System Lab
