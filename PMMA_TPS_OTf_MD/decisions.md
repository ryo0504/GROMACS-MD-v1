# 今回の前提・採用条件・未確定事項

## 2026-10-28の依頼で確定した内容

- 今回はMacでの入力準備まで。実分子QM/RESP・力場検証・MDは研究室のLinuxサーバーで行う。
- 最終的な目標は、検証済みのMD入力一式まで完成させるA案。
- Gaussian 16が使用可能。ESPはMK格子を採用するという利用者回答がある。
- g16に使えるCPU・メモリとPBS例は不明。明日確認する。
- 保存先は利用者が作成した`PMMA_TPS_OTf_MD/`。v2仕様書や既存エタノール教材を変更しない。

## v2との関係

| 項目 | v2 | 今回の準備 |
|---|---|---|
| QM backend | PySCFが第一候補 | Gaussian 16へ変更 |
| QM法・基底 | 気相RHF/6-31G* | 維持 |
| ESP | 独自Fibonacci shell | **Gaussian MKへ変更（利用者指定）** |
| 電荷fit | RESP、内部/末端template、独立評価 | 維持 |
| 部分電荷のscale | 1.0 | 維持 |
| 分子・個数・箱 | PMMA19 100/TPS83/OTf83、10 nm角 | 維持 |
| 高温MD | 1400 K、100→10→1 bar、0.5/0.5/5/20 ns | 維持 |
| 冷却・室温 | 主経路77.2 ns、室温NVT追加で97.4 ns | 維持 |
| 保存構成 | repo内の複数領域 | 指定フォルダ内へ集約して移送しやすくした |
| 設定ファイル | YAML例 | 標準Pythonで読むJSON。内容・役割は同じ |

MK点列はv2の点列と同一ではないため、二つを同じfitに混ぜない。GaussianのMKで表示されるESP fit電荷を、そのままRESP電荷と呼ばない。

## 明日確認する項目

1. PBS queue、許可されたCPU/メモリ/walltime、Gaussianの使用ルール。
2. g16の環境有効化、`formchk`、scratchのパスと容量。
3. Pythonの実体、NumPy/RDKit/ParmEd/pyedr、AmberTools/ACPYPE/GROMACSの実体。
4. g16のMK ESPを採用AmberToolsの`espgen`が正しく読めるか。
5. 最初のOTfジョブの所要時間・メモリ・収束。その実測を基に残りの資源を決める。

## 未実施

実分子Gaussian、ESP抽出、RESP、電荷転写の物理的評価、GAFF2、変換エネルギー比較、packing、grompp/MD、PBS投入は未実施。人工ESPでの拘束テストは実分子の電荷計算とは分離する。
