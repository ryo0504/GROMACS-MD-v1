# Python environment と host binary

`environment.yml` は科学Pythonのsource of truth。Python3.11を選び、platform固有buildはexplicit lockへ記録します。
AmberTools24とNumPy1.26系を共通の安定した基盤とします。最新版へ無条件に追従しません。
MACEのmatscipyはNumPy1系と互換な1.1.1を選択します。
`setup_micromamba.sh` は対応platformの `env/locks/<platform>.lock` を優先し、なければYAMLをsolveします。
lockはConda packageだけです。MACE/PyTorch/pipの実際のversionは追加inventoryにも記録します。
`<platform>-pip.txt` はConda管理外の補助packageだけをversion固定します。対応するConda lockから環境を作成した後、
同じbackendの再構築なら `MACE_PIP_LOCK=env/locks/osx-arm64-pip.txt ./scripts/setup/install_mace.sh` が使えます。
Linuxでは同じlinux-64でもCPU/CUDA backendが違うため、backendを確認してから明示的にlockを選択してください。
portable YAMLへホストのbinary buildを大量に固定しません。

MacのPyTorchは同じConda runtimeを使うconda-forge版で導入します。pip wheelと別のOpenMP runtimeを混在させることを避けます。
Linux NVIDIAの場合は `nvidia-smi` とserver policyを確認し、公式 [PyTorch selector](https://pytorch.org/get-started/locally/) に従って
`TORCH_VERSION` / `TORCH_INDEX_URL` を指定します。nvccとPyTorch wheelのruntimeは同一である必要はありませんがdriver互換性が必要です。
CPU-only LinuxはCPU wheelを利用します。既存PyTorchがあれば維持し、NVIDIA hostではcuda利用可能性を検証します。

```bash
TORCH_VERSION=<REVIEWED_VERSION> TORCH_INDEX_URL=<OFFICIAL_INDEX> ./scripts/setup/install_mace.sh
./scripts/setup/export_environment.sh
```

GROMACSはenvironment.ymlに含めません。既存command/moduleが優先です。Macローカルbuildはthread-MPI/CPU版でOpenMPはOFF、Apple clangを使用します。
Linux CUDA/external MPIはprobe後に決めます。
同じ科学入力でも異なるbinary/GPU/precisionでtrajectoryがbitwise一致することは保証しません。

## このPCの実環境

2026-10-05に環境を同期対象外へ再構築しました。実パス、カーネル、確認結果は [current_pc_setup](../docs/current_pc_setup.md) を参照してください。scripts/setup/common.sh は `.local/host.env` がある場合に読み込みます。明示されたMAMBA_ROOT_PREFIX/GMX_BINを優先します。

## 無料のローカルQM / RESP教材（任意追加）

エタノールNotebookの第5節はPySCF 2.14.0とgeomeTRIC 1.1.1を使用します。既存`gmx-md`へ追加する場合は、repository rootで以下を実行し、Notebookのカーネルを再起動します。

```bash
source scripts/setup/common.sh
"$MM" run -n gmx-md python -m pip install --prefer-binary -r env/requirements-qm.txt
"$MM" run -n gmx-md python -m pip check
"$MM" run -n gmx-md python -c 'import pyscf, geometric; print(pyscf.__version__, geometric.__version__)'
```

この追加は`env/requirements-qm.txt`に固定し、従来のConda lockを変更しません。古いlockで環境を再構築した場合は、上の追加操作も必要です。PySCF wheelは計算コードを含むため、純Pythonのみの実装ではありません。[PySCF公式導入手順](https://pyscf.org/user/install.html)も参照してください。
