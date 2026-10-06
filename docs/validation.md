# v1 検証記録

更新：2026-10-05にこのPCの環境を再構築し、GROMACSとgmx-mdカーネルを使用可能にしました。最新の状態は [環境構築記録](current_pc_setup.md) を参照してください。以下は初回実装時の検証記録です。

検証日: 2026-10-04（Asia/Tokyo）。元の仕様書・指示書と、作業中に追加されたEthanol benchmark仕様書を保持しました。
今回の実装範囲はGROMACS MD Workbench v1です。別仕様の27分子NVT benchmarkは実装・計算していません。

## 実機

| 項目 | 状態 |
|---|---|
| OS / CPU | macOS 26.6.2 / arm64 |
| System Python | /usr/bin/python3 3.9.6、変更なし |
| Workbench Python | Micromamba gmx-md: 3.11.16 |
| Micromamba | 2.9.0、.local/bin |
| GROMACS | 2025.4、.local/gromacs/bin/gmx |
| GROMACS build | CPU / mixed precision / ARM_NEON_ASIMD / thread-MPI / OpenMP OFF / GPU OFF |
| AmberTools | 24.8 |
| ACPYPE | 2023.10.27 |
| ASE / MACE | 3.29.0 / 0.3.16 |
| PyTorch | conda-forge 2.8.0、CPU/MPS available、CUDA unavailable |
| NumPy / MDTraj | 1.26.4 / 1.10.3 |
| Homebrew | 未導入 |
| Jupyter kernel | Python (gmx-md)をuser kernelとして登録 |

## 結果

- pytest: **39 passed**（12.91秒、skipなし）。
- pip check: **No broken requirements found.**
- RDKit構造生成→ASE変換→SDF/XYZ保存: 成功。
- AmberTools antechamber→parmchk2→tleap→ACPYPE: 実commandで成功。
- GROMACS grompp/mdrun: EM→NVT→NPT→productionの全4段階成功、trajectory生成確認。
- Run入力改変、同じstage出力への再実行、missing RESP sourceの拒否: 確認。
- 4 Notebook: nbformat検証とgmx-md kernelでの実行成功。
- doctor: 全16 Python moduleをimport、全AmberTools executable/ACPYPE/GROMACSを検出。
- setup再実行: 既存環境を保持しpackage/kernel登録を再確認。
- shell syntax / Python compile: 成功。
- Git ignore: 大容量出力/環境を除外、小さい科学入力は管理可能。
- Linux-64環境: Linux 5.4/glibc2.28/CPU-onlyを**仮想条件として指定したdry-run solve**成功。server実機での導入・GPU/MPI/SLURM実行は未検証。

integration testのchargeは合成テストデータです。実際のQM/RESP結果ではなく、科学的なethanol simulationの検証を意味しません。
MACEはcalculator importまで確認しました。利用者の科学的modelが未指定のためenergy/force/optimizationの実model検証は行っていません。
Gaussian ESP→RESP interfaceは用意しましたが、実QM/ESPとrestraint inputがないため実データfitは未実行です。

## 依存関係・ビルドの判断

最新AmberTools26で付属ツールのNumPy要件に不整合があったため、AmberTools24とNumPy1.26、対応MDTraj/Matscipyを選択しました。
Macでpip PyTorchとConda OpenMP runtimeの衝突を実際に確認し、Conda版PyTorchへ統一して解消しました。
unsafeなKMP_DUPLICATE_LIB_OKは設定していません。
GROMACS付属FFTWのAutotoolsがdecomposed UnicodeのVPATHで失敗したため、ASCIIのtemporary directoryでbuildし、リポジトリ内prefixへinstallしました。

GROMACS source: [official archive](https://ftp.gromacs.org/gromacs/gromacs-2025.4.tar.gz)

```text
SHA256 ca17720b4a260eb73649211e9f6a940ee7543452129844213c3accb0a927a5c3
```

Conda build lockとpip補助package versionは `env/locks/osx-arm64*` に保存しています。
gitはmain branchで初期化済み、commit/remoteは未設定です。commitがないrunではmetadataのcommitはnullです。
