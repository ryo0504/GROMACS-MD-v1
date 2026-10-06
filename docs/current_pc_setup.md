# このPCの環境構築記録

構築日：2026-10-05（JST）。macOS26.6.2 / Apple Silicon arm64。

## 使用するカーネル

**Python (GROMACS MD Workbench)**、カーネルID **gmx-md**。

VS CodeではNotebook右上の「カーネルの選択」からこの名前を選びます。表示が更新されない場合はNotebookを開き直し、「別のカーネルを選択」→「Jupyterカーネル」から選びます。Python interpreterを直接指定する場合：

```text
$HOME/.local/share/mdworkbench/GROMACS-MD-v1/micromamba/envs/gmx-md/bin/python
```

カーネルの実登録先：

```text
$HOME/Library/Jupyter/kernels/gmx-md/kernel.json
```

このカーネルはPythonだけでなく、AmberToolsのPATH/AMBERHOMEと外部GROMACSのパスも保持します。VS CodeをDockから開いた場合でも外部commandを検出できるようにしました。再登録は `scripts/setup/register_kernel.py` をgmx-md環境内で実行します。

## すぐに使う

repository rootで：

```bash
./scripts/setup/start_jupyter.sh
```

これはgmx-md環境のJupyterLabでエタノール教材を開きます。terminalでも環境を有効にする場合：

```bash
source scripts/setup/common.sh
eval "$("$MM" shell hook -s zsh)"
micromamba activate gmx-md
python --version
gmx --version
make doctor
```

activateせずにcommandを実行する場合：

```bash
source scripts/setup/common.sh
"$MM" run -n gmx-md python -m mdworkbench.ethanol_benchmark build --output runs/ethanol27_initial_02
```

build先を新しい名前にしたらNotebookのBUILDも合わせます。既存runを上書きしません。

## 構築内容

| 項目 | 確認結果 |
|---|---|
| Micromamba | 2.9.0 / osx-arm64 |
| Python | 3.11.16 |
| NumPy / ASE | 1.26.4 / 3.29.0 |
| AmberTools / ACPYPE | 24.8 / 2023.10.27 |
| MACE / PyTorch | 0.3.16 / 2.8.0 |
| PyTorch device | CPU使用可、MPS使用可、CUDAなし |
| GROMACS | 2025.4 / mixed precision / ARM_NEON_ASIMD |
| GROMACS parallelism | thread-MPI、OpenMP OFF、GPU OFF |
| Workbench | editable install済み |
| JupyterLab / ipykernel | 導入・カーネル登録済み |

環境本体はGoogle Driveの同期対象外に配置しました：

```text
$HOME/.local/share/mdworkbench/GROMACS-MD-v1/
├── bin/micromamba
├── micromamba/envs/gmx-md/
└── gromacs/
    ├── bin/gmx
    ├── lib/
    └── share/gromacs/top/
```

repositoryの `.local/host.env` はこのPCの参照先を持ち、`scripts/setup/common.sh` が読み込みます。このファイルはGit対象外です。参照先のMicromamba binaryが存在しないPCでは設定を適用しません。terminalで明示したMAMBA_ROOT_PREFIX/GMX_BINは優先します。

元のrepository内 `.micromamba` / `.local/gromacs` は保持していますが、現在の実行には上記のPC専用環境を使います。system Python、Homebrew、shell rc、CUDA driverは変更していません。

## 修復と再現性

同期された環境にある共有ライブラリへのsymlinkが欠落していました。Python環境は既存 `env/locks/osx-arm64.lock` と `osx-arm64-pip.txt` から、このPCのMicromambaで再構成しました。共有環境のnative libraryはconda-forge版に揃え、別のOpenMP runtimeを重複ロードする回避変数は使っていません。

GROMACSの既存arm64 binaryをPCローカルへ配置し、`libgromacs.10.dylib` と `libmuparser.2.dylib` の参照を実体へ戻しました。force-field dataは同versionの公式sourceから配置し、Data prefixが新しいローカルprefixになることを確認しました。再ビルドは不要でした。

公式sourceのSHA256：

```text
ca17720b4a260eb73649211e9f6a940ee7543452129844213c3accb0a927a5c3
```

再現用Conda lockとpip inventoryは `env/locks/osx-arm64*` に保存しました。別のMacやLinuxでは、そのマシンで環境を再構成します。

## 確認結果

- pip check：No broken requirements found。
- doctor：全16 Python moduleをimport、AmberToolsの5 command、ACPYPE、GROMACSを検出。
- pytest：**46 passed、skipなし**。AmberTools→ACPYPEとGROMACS EM→NVT→NPT→productionも実commandで成功。
- 全5冊のNotebookを、登録したgmx-mdカーネルで実行。コードセルのerrorなし。
- shell構文、Notebook schema、kernel.jsonの実Pythonパスを確認。

ログ・実行済みNotebookは `runs/environment-setup-20261005/` に保存し、Git対象外です。

integration testの電荷は合成fixtureであり、科学的RESP結果ではありません。エタノール27分子教材のRESP/QM入力や力場条件はまだ未準備です。この環境構築はその値を変更していません。MACEはcalculator importまで確認し、未指定modelの計算は開始していません。

公式資料：[Micromamba](https://mamba.readthedocs.io/en/stable/user_guide/micromamba.html)、[GROMACS 2025.4 installation](https://manual.gromacs.org/documentation/2025.4/install-guide/index.html)、[GROMACS source archive](https://ftp.gromacs.org/gromacs/gromacs-2025.4.tar.gz)。
