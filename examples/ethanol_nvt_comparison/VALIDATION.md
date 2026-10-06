# 検証記録（2026-10-05 JST）

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
