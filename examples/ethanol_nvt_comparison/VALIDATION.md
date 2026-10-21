# 検証記録

## 2026-10-21：PySCF → RESP → GAFF2 → 実NVT

同じ`gmx-md` / `Python (GROMACS MD Workbench)`カーネルで、1冊の`ethanol_comparison.ipynb`を更新し、全18コードセルを実行しました。再実行でも全セルが成功し、GROMACSのcommand履歴が増えず、完了済み計算をhash検証後に再利用できることを確認しました。

| 検証項目 | 実際の結果 |
|---|---|
| PySCF / geomeTRIC | 2.14.0 / 1.1.1。H2 RHF/STO-3Gの収束とimport確認済み |
| 気相エタノール | RHF/6-31G*、構造最適化・最終SCFとも収束、−154.0743658522 Hartree |
| ESP | 4層の独自露出Fibonacci表面格子625点。核電荷と電子密度から計算 |
| 2段階RESP | 両stage収束、C結合Hの群内等価性とO/OHのHの固定を検証 |
| RESP電荷総和 | 約+0.000001 e（6桁丸め）。AmberToolsの中性化でOに−0.000001 eの差、topology総和は数値誤差内で0 |
| ESP RMS誤差 | 0.00175338 Hartree/e（fitに使用した625点で評価） |
| 力場 | GAFF2 → Amber topology → ACPYPEを実commandで実行 |
| 原子mapping | 座標・元素・結合・電荷を照合、topology→ASEは0〜8の順序を保持 |
| 相互作用 | 8 bonds / 13 angles / 12 proper dihedrals / 12 explicit 1–4 pairs、nrexcl=3 |
| 組合せ・1–4 | combination rule 2、LJ係数0.5、Coulomb係数1/1.2 |
| 非結合条件 | PME、LJ Cut-off、rlist/rcoulomb/rvdw=0.60 nm、EnerPres。採用理由をconfigへ記録 |
| ゼロステップ評価 | safe=true、Potential −1228.199707 kJ/mol、最大force 2008.6782 kJ/mol/nm、最短分子間距離1.7446977 Å |
| MD | 600 steps、0〜0.3 ps、243原子、固定箱、有限座標、XTC 121 frames |
| 解析 | step205〜600の80frames。温度・PotentialとOO/OH RDFの図を表示・目視確認 |
| 初期条件保存 | 元のinitial.npz / initial.extxyzのhash一致。GROへの座標・速度転記も照合 |
| 検査 | pytest **52 passed**、skipなし。pip check正常、Notebook schema正常、git diff --check正常 |

**限定的な警告許可**：GROMACS 2025.4が意図した`comm-mode=none`に対して出す既知のCOM警告は、2026-10-21にユーザーが採用した限定例外に従って扱いました。通常のgromppを先に実行し、その警告1件だけが停止原因であると内容照合した場合だけ`-maxwarn 1`を付けて再実行。再実行後も警告内容を照合し、両command/logと理由をmetadataへ保存しました。別の警告やエラーでは停止することをテストしています。汎用Workbenchの挙動は変更していません。

追加テストでは、正規化s型Gaussianに電子1個を置いた解析的な静電ポテンシャルとPySCF積分の結果を比較し、核・電子の符号とbohr単位を検査しました。Amber ESPの固定幅入出力、表面格子の露出条件、hash改変時の停止も検証しました。

保存先：`runs/tutorial-pyscf-20261021/`（実行済みNotebook・図・検証資料）、`runs/ethanol27_pyscf/`（QM/RESP/config）、`runs/ethanol27_pyscf_parameters/`（力場）、`runs/ethanol27_pyscf_nvt/`（MD/解析）。元の未設定テンプレート`systems/ethanol_27/system.yaml`は保持しました。

**評価の限界**：独自ESP格子はGaussianのMK格子と同一ではありません。単一気相配座、振動解析なし、小さいcutoff、300 fsの操作教材です。温度が初期に大きく変動し、Potentialも後半まで変化しており、平衡・RDF収束・実験物性との一致を主張しません。MACEの比較元は引き続き未準備で、MACE計算・比較は実施していません。

## 2026-10-05：初回実装時の記録（以下は当時の状態）

更新：2026-10-05にこのPCの環境を再構築し、GROMACSとgmx-mdカーネルを使用可能にしました。最新の状態は [環境構築記録](../../docs/current_pc_setup.md) を参照してください。以下は初回実装時の検証記録です。

## 作成状況

作成前はこの仕様に対応するsystems/ethanol_27とexamples/ethanol_nvt_comparisonはなく、汎用Workbenchの4冊のみでした。今回専用教材を新規作成しました。

## 実行して確認した範囲

- repository全体のpytest：**38 passed, 8 skipped**。未導入の外部softwareや任意packageに依存する項目はskip。
- Notebookの全10コードセルを実行。実MD出力がない場合は未実行の表示で止まり、架空の結果を作らない。
- ASE3.29.0、NumPy1.26.4、Python3.11.14で構造生成。
- 27分子、243原子、1分子46.069 amu、volume2628.072 Å³、density0.785930943 g/cm³。
- independent seed42 placement/velocities、system COM≈0、726DoFの300 K、729DoFの298.765432 K。
- 初期構造、速度分布、初期RDFのPNGを確認。軸・単位・legendと箱表示に問題なし。
- Notebook schema validation、Python compile、追加ファイルの空白検査。
- RDFの分子内pair除外、MIC、r_max、45bins、26/V、step205〜600の80枚をテスト。
- 明示的な科学入力がない場合にprepareが停止することをCLIで確認。
- **合成のゼロ電荷test fixtureのみ**で27分子topology展開、GROの高精度座標・速度保持、input hash改変検出を検証。これはRESP/GAFF科学計算の成功を意味しない。fixtureはpytestの一時directoryだけに生成。
- ParmEdが不足する1–4pairを警告した場合、prepareは停止する。
- Temperature/PotentialのXVG列順をlegendから解釈するテスト。

## 検証環境と成果物

既存 `.micromamba/envs/gmx-md` Pythonは当初実行権限がなく、実行権限を付けて起動を試した後もimportに長時間を要しました。既存GROMACSも実行権限を付けて起動を試しましたが、`libgromacs.10.dylib` の参照が見つからず起動できませんでした。既存environmentの再構築やライブラリ書換えは行っていません。

検証は一時的な `/tmp/ethanol27-py311` 環境を使用。Micromamba管理の既存環境へのpackage追加はしていません。`env/ethanol27-study.yml` は別の学習用Micromamba環境を作る定義ですが、そのsolve/installは未検証です。

出力はrepository rootの `runs/ethanol27_initial/` と `runs/ethanol27_validation/`。実行済みNotebookは `runs/ethanol27_validation/ethanol_comparison_executed.ipynb`。実package inventoryは `runs/ethanol27_initial/validation_environment.txt`。これらはGit対象外です。

## 未実行・未確定

- QM/ESP/RESP、科学的GAFF/GAFF2単分子topologyの作成。
- GAFF/GAFF2選択、cutoff/長距離処理/dispersion correction/force thresholdの科学判断。
- 実GROMACSのgrompp、ゼロstep preflight、600-step NVT、EDR/XTC解析。
- MACE比較元はないため、元の初期座標・速度との一致とMACE/GROMACS比較。
- MACE0.3.16実行環境の再現。今回はMACE実行不要でありversion未導入。

仕様書の全完成条件を達成した、あるいはMDを実行済みとは報告しません。教材の前半と実行/解析interfaceを作成・検証済みです。
