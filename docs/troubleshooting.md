# Troubleshooting

| 症状 | 対応 |
|---|---|
| micromamba command not found | `.local/bin`をPATHへ追加、setupを実行 |
| Python3.9が使われる | shell hookを評価しgmx-mdをactivate、which pythonを確認 |
| GROMACS未設定 | 既存module/GMXRCを設定、Macはinstall_gromacs_mac.sh |
| Mac compiler不足 | 利用者が `xcode-select --install` を実行 |
| RESP source不足 | 実際のQM/ESP/fitted chargesを保存してconfigを指定 |
| missing include | topologyとincludeを保存、prepare-runへ--include-dir |
| grompp warning/error | stage/grompp.logとMDP/systemを修正。-maxwarnで隠さない |
| already existing output | 新しいrun/build/parameter directory名を選択 |
| stale .running lock | running processを確認した後、利用者が空lock directoryをrmdir |
| CUDA unavailable | driver、torch backend、SLURM GPU allocationを確認 |
| OpenMP duplicate runtime | conda-forge PyTorch/runtimeに統一。KMP_DUPLICATE_LIB_OKは使わない |
| Notebookのmodule不足 | kernelをPython (gmx-md)に変更 |

RESP atom countが合っていても順序やconformationが合うとは限りません。実際の対応を確認してください。
テストでdependencyが未導入ならskipしますが、導入済みdependencyが壊れている場合は失敗として修正してください。
