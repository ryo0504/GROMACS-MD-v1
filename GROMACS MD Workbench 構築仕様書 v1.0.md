# GROMACS MD Workbench 構築仕様書 v1.0

## 1. 目的

GROMACSを中心とした分子動力学計算環境を構築する。

開発・入力作成・小規模テストはMacで行い、本番計算は将来的にLinux計算サーバーで実行する。

GitHubを介して、

- ソースコード
- Notebook
- 計算条件
- 力場パラメータ
- 環境定義
- 実行スクリプト
- 解析コード

を共有し、異なるマシン上でも再現可能なMDワークフローを構築する。

仮想環境そのものはGit管理せず、「環境を再構築するための定義」をGit管理する。


# 2. 基本方針

## 2.1 開発機と本番機の役割

### Mac

用途：

- VS Codeによる開発
- Jupyterによる探索
- RDKitによる分子構造生成
- ASEによる構造操作
- MACEによる補助的な構造最適化・energy/force評価
- AmberToolsによるGAFF系パラメータ作成
- GROMACS input生成
- 小規模GROMACSテスト
- MD trajectory解析
- 可視化

Macでは計算性能を最大化することを目的とせず、

「本番計算前にworkflowが正しく動くことを確認する環境」

と位置付ける。


### Linux計算サーバー

用途：

- GROMACS本計算
- GPU/CPUを利用した長時間MD
- 必要に応じてMPIによる並列計算
- SLURM等によるジョブ投入
- trajectory保存

Linuxではハードウェアに合わせてGROMACSを最適化する。


# 3. 環境管理方針

## 3.1 Python環境

Micromambaを第一選択とする。

Python関連ツールは原則として1つのMicromamba環境

`gmx-md`

にまとめる。

主な対象：

- Python
- JupyterLab
- ipykernel
- NumPy
- SciPy
- pandas
- matplotlib
- RDKit
- ASE
- MDAnalysis
- MDTraj
- ParmEd
- AmberTools
- ACPYPE
- PyYAML
- pytest

MACEについてはPyTorch依存性がプラットフォーム・GPUによって変化するため、同じ `gmx-md` 環境へプラットフォーム依存インストールする。


## 3.2 GROMACS

GROMACS本体はMicromamba環境から分離する。

理由：

- CPUアーキテクチャ依存
- SIMD依存
- GPU依存
- CUDA依存
- MPI依存
- HPC環境依存

MacとLinuxで同一バイナリ環境を再現することは目標としない。

代わりに、

- GROMACS version
- build information
- GPU backend
- MPI availability
- executable path

を記録する。


# 4. 再現性の考え方

完全なバイナリ環境一致ではなく、

「同じ科学的条件を再構築できること」

を再現性の目標とする。

Git管理対象：

- Python source
- Notebook
- YAML/TOML設定
- MDP
- 分子構造
- RESP charge
- GAFF/GAFF2 parameter
- top/itp
- 実行script
- SLURM template
- environment定義
- README
- validation code

Git管理対象外：

- trajectory
- checkpoint
- 大容量energy file
- temporary file
- cache
- 仮想環境
- MACE model cache
- Jupyter checkpoint
- GROMACS backup
- 大規模log
- GPU cache


# 5. Repository構造

以下を基本構造とする。

```text
gromacs-md-workbench/
│
├── README.md
├── pyproject.toml
├── Makefile
├── .gitignore
├── .gitattributes
│
├── env/
│   ├── environment.yml
│   ├── requirements-mace.txt
│   ├── README.md
│   └── locks/
│       └── .gitkeep
│
├── config/
│   ├── defaults.yaml
│   ├── example_system.yaml
│   └── platform.example.yaml
│
├── src/
│   └── mdworkbench/
│       ├── __init__.py
│       ├── chemistry/
│       │   ├── build.py
│       │   └── convert.py
│       │
│       ├── parameterization/
│       │   ├── amber.py
│       │   ├── resp.py
│       │   └── validation.py
│       │
│       ├── gromacs/
│       │   ├── prepare.py
│       │   ├── runner.py
│       │   └── validation.py
│       │
│       ├── analysis/
│       │   ├── trajectory.py
│       │   └── rdf.py
│       │
│       ├── provenance/
│       │   └── metadata.py
│       │
│       └── utils/
│           └── subprocess.py
│
├── notebooks/
│   ├── 00_environment_check.ipynb
│   ├── 01_build_molecule.ipynb
│   ├── 02_parameterization_review.ipynb
│   └── 03_trajectory_analysis.ipynb
│
├── mdp/
│   ├── em.mdp
│   ├── nvt.mdp
│   ├── npt.mdp
│   └── md.mdp
│
├── scripts/
│   ├── setup/
│   │   ├── doctor.sh
│   │   ├── setup_micromamba.sh
│   │   ├── install_mace.sh
│   │   ├── detect_gromacs.sh
│   │   └── probe_server.sh
│   │
│   ├── prepare/
│   │   ├── build_system.py
│   │   └── parameterize.py
│   │
│   └── run/
│       ├── em.sh
│       ├── nvt.sh
│       ├── npt.sh
│       └── production.sh
│
├── cluster/
│   ├── README.md
│   ├── slurm/
│   │   ├── cpu.slurm.template
│   │   └── gpu.slurm.template
│   └── detected/
│       └── .gitkeep
│
├── systems/
│   └── example/
│       ├── README.md
│       ├── system.yaml
│       ├── structures/
│       └── parameters/
│
├── runs/
│   └── .gitkeep
│
├── tests/
│   ├── test_imports.py
│   ├── test_structure.py
│   ├── test_parameterization.py
│   └── test_gromacs_detection.py
│
└── docs/
    ├── architecture.md
    ├── workflow.md
    ├── parameterization.md
    ├── server_setup.md
    └── troubleshooting.md
```


# 6. Notebook / Python / Shell の役割

以下の責務分離を厳守する。

## Notebook

用途：

- 探索
- 分子確認
- 可視化
- trajectory解析
- グラフ作成
- interactive analysis

Notebook内に本番workflowのロジックを大量に書かない。


## Python script / src

用途：

- 再利用可能な処理
- 分子生成
- file conversion
- parameterization
- validation
- GROMACS input作成
- metadata生成
- trajectory解析

科学的処理の本体はこちらに置く。


## Shell / SLURM

用途：

- GROMACS実行
- HPC job submission
- environment activation
- resource allocation

本番計算はNotebookから直接実行しない。


# 7. 分子構造生成

初期対象：

- 有機低分子
- 簡単なオリゴマー
- 簡単な高分子モデル
- 多成分有機系

無機結晶への対応はv1では不要。


## RDKit

主用途：

- SMILES読み込み
- H付加
- 3D conformer生成
- 基本的な構造操作
- SDF/MOL出力


## ASE

主用途：

- xyz操作
- geometry manipulation
- calculator interface
- format conversion


## MACE

主用途：

- geometry optimization補助
- energy評価
- force評価

MACEはGROMACSの力場として使用しない。

v1では

`MACE MD`

と

`GROMACS classical MD`

を明確に分離する。


# 8. GAFF / RESP parameterization

主な力場候補：

- GAFF
- GAFF2

charge model：

- RESP

を想定する。


## 8.1 重要原則

RESP chargeは科学的入力条件である。

以下をAIが勝手に決定しない。

- QM method
- basis set
- molecular conformation
- total charge
- multiplicity
- ESP計算方法
- RESP fitting条件

これらは `system.yaml` に明示する。


## 8.2 system.yaml例

```yaml
system:
  name: ethanol
  net_charge: 0
  multiplicity: 1

force_field:
  family: gaff2
  charge_model: resp

resp:
  qm_backend: external
  method: HF
  basis: 6-31G*
  source_file: null

mace:
  use_geometry_optimization: false

gromacs:
  temperature_K: 298.15
  pressure_bar: 1.0
```


## 8.3 RESP workflow

想定workflow：

```text
SMILES / SDF / MOL2
        │
        ▼
RDKit structure generation
        │
        ▼
optional geometry optimization
        │
        ▼
QM ESP calculation
        │
        ▼
RESP fitting
        │
        ▼
GAFF / GAFF2 atom typing
        │
        ▼
parmchk2
        │
        ▼
tleap
        │
        ▼
AMBER topology
        │
        ▼
ACPYPE / validated conversion
        │
        ▼
GROMACS
.top / .itp / .gro
```

v1ではQM calculation backendを必須実装しなくてもよい。

代わりに、

「外部で生成されたQM output / ESPを受け取って続きを実行できる」

構造を優先する。

RESP sourceが存在しない場合に、勝手にAM1-BCC等へfallbackしてはならない。

明示的なerrorを出すこと。


# 9. Parameter validation

GROMACS実行前に以下を検証する。

- net charge
- atom count
- residue count
- atom type
- missing parameter
- topology include
- coordinate/topology consistency
- molecule count
- periodic box
- GROMACS grompp

parameterization時には使用した、

- AmberTools version
- ACPYPE version
- GAFF/GAFF2
- charge model
- RESP条件
- source structure hash

をmetadataへ保存する。


# 10. GROMACS workflow

基本workflow：

```text
structure
   │
   ▼
parameterization
   │
   ▼
system construction
   │
   ▼
grompp
   │
   ▼
Energy Minimization
   │
   ▼
NVT
   │
   ▼
NPT
   │
   ▼
Production MD
   │
   ▼
Trajectory analysis
```

MDPは

- `em.mdp`
- `nvt.mdp`
- `npt.mdp`
- `md.mdp`

に分離する。


# 11. run directory

計算ごとに独立directoryを作る。

例：

```text
runs/
└── ethanol_20261004_180000/
    ├── metadata.json
    ├── input/
    ├── em/
    ├── nvt/
    ├── npt/
    └── production/
```

`runs/` はGit管理しない。


# 12. Provenance

各run開始時に `metadata.json` を自動生成する。

最低限以下を記録する。

```text
timestamp
hostname
OS
CPU architecture
GPU
git commit hash
git dirty status
Python version
Micromamba environment
installed Python packages
GROMACS version
GROMACS executable
GROMACS build information
AmberTools version
ACPYPE version
MACE version
PyTorch version
CUDA availability
input file hashes
MDP hashes
topology hashes
random seed
executed command
```

これにより、

「いつ、どのコード、どの環境、どのinputで計算したか」

を追跡可能にする。


# 13. .gitignore

最低限以下を除外する。

```gitignore
# Python
__pycache__/
*.py[cod]
.pytest_cache/
.ruff_cache/

# Virtual environments
.venv/
venv/
.micromamba/
micromamba/

# Jupyter
.ipynb_checkpoints/

# Local machine configuration
.local/
.env
config/platform.local.yaml

# GROMACS outputs
runs/
*.xtc
*.trr
*.tng
*.edr
*.cpt
*.log
*.tpr
*.xvg

# GROMACS backups
#*
*#

# Large temporary structures
*.tmp
tmp/
temp/

# MACE / PyTorch
*.model
*.pt
*.pth
.cache/
torch_cache/

# macOS
.DS_Store
```

ただし、

- `.top`
- `.itp`
- `.mdp`
- 最終parameter `.mol2`
- RESP charge情報
- 小さい `.gro`
- reproducibilityに必要なparameter file

は原則Git管理可能とする。


# 14. Micromamba

`env/environment.yml` を環境定義のsource of truthとする。

例：

```yaml
name: gmx-md

channels:
  - conda-forge

dependencies:
  - python=3.11
  - pip
  - jupyterlab
  - ipykernel
  - numpy
  - scipy
  - pandas
  - matplotlib
  - rdkit
  - ase
  - mdanalysis
  - mdtraj
  - parmed
  - ambertools
  - acpype
  - pyyaml
  - pytest
```

厳密なbuild番号を `environment.yml` に大量に固定しない。

実際に使用したプラットフォームについてはlock/spec生成を可能にする。

例：

```text
env/locks/osx-arm64-*.lock
env/locks/osx-64-*.lock
env/locks/linux-64-*.lock
```

利用可能なlockがあれば優先して使用する。

存在しないプラットフォームでは `environment.yml` からsolveし、新しいlockを生成できるようにする。


# 15. MACE installation

MACEは `gmx-md` Micromamba environment内部へインストールする。

ただしPyTorch backendを固定しない。

`scripts/setup/install_mace.sh`

は以下を検出する。

```text
OS
architecture
CUDA
NVIDIA GPU
PyTorch availability
```

Mac：

- CPU/MPS対応PyTorch
- mace-torch

Linux NVIDIA：

- server CUDA環境と互換なPyTorch
- mace-torch

を選択する。

CUDA versionを推測だけで決めない。


# 16. GROMACS detection

`scripts/setup/detect_gromacs.sh`

を作成する。

確認項目：

```bash
which gmx
gmx --version
```

必要に応じて：

```bash
which gmx_mpi
gmx_mpi --version
```

GROMACSが存在しない場合、

Python環境構築全体を失敗扱いにはしない。

代わりに、

「GROMACS未設定」

と明確に報告する。


# 17. Mac GROMACS

Macではテスト実行を目的とする。

ローカルGROMACSはsystem側に置く。

Micromamba environmentへGROMACS本体を含めない。

既存のGROMACSを最優先する。

存在しない場合はMacのCPU architectureと利用可能なpackage managerを検出し、インストール方法を提示・スクリプト化する。

GROMACSの性能最適化より、workflow検証を優先する。


# 18. Linux server adaptation

本番server仕様は現時点で固定しない。

その代わり、

`scripts/setup/probe_server.sh`

を作成する。


## probe項目

最低限以下を調査する。

```text
Linux distribution
kernel
architecture
CPU model
physical cores
logical cores
RAM
GPU
nvidia-smi
NVIDIA driver
CUDA
nvcc
gcc
clang
cmake
MPI
mpirun
OpenMPI
SLURM
sbatch
srun
module command
existing GROMACS
filesystem
home quota if obtainable
scratch filesystem if obtainable
```

結果を

```text
cluster/detected/platform_report.txt
```

へ保存する。


# 19. 本番serverでAIが行う判断

platform reportを読んだCodex/Copilotは、

1. 既存moduleとしてGROMACSが提供されていないか確認
2. 既存GROMACSが十分なら再利用
3. NVIDIA GPUならCUDA backendを検討
4. single nodeならthread-MPIも候補
5. multi-node利用ならexternal MPIを検討
6. SLURM環境ならSLURM job scriptを生成
7. sudo権限を前提としない
8. cluster policyを尊重
9. system-wide packageを勝手に変更しない

こと。


# 20. SLURM

v1ではtemplateのみ作成する。

```text
cluster/slurm/cpu.slurm.template
cluster/slurm/gpu.slurm.template
```

実際の

- partition
- account
- number of CPUs
- GPU type
- memory
- walltime

はserver環境確認後に設定する。


# 21. Doctor command

以下を一括確認できるようにする。

```bash
./scripts/setup/doctor.sh
```

出力例：

```text
[OK] micromamba
[OK] Python
[OK] RDKit
[OK] ASE
[OK] AmberTools
[OK] ACPYPE
[OK] MDAnalysis
[OK] Jupyter
[OK] MACE
[OK] PyTorch
[INFO] CUDA unavailable
[OK] GROMACS 20xx.x
```

エラーの場合は、

「何が足りないか」
「どう修正するか」

を表示する。


# 22. Jupyter kernel

Micromamba環境をJupyter kernelとして登録する。

kernel name：

```text
gmx-md
```

Notebookから必ずこのkernelを選択可能にする。


# 23. Test

pytestを導入する。

最低限確認：

- RDKit import
- ASE import
- MDAnalysis import
- ParmEd import
- AmberTools executable detection
- ACPYPE detection
- MACE import
- PyTorch import
- GROMACS detection
- structure conversion
- configuration parser
- metadata generation

外部GROMACSがない場合はGROMACS integration testをskip可能にする。


# 24. Notebook

## 00_environment_check.ipynb

以下を表示：

- Python
- RDKit
- ASE
- MACE
- torch device
- AmberTools
- ACPYPE
- MDAnalysis
- GROMACS


## 01_build_molecule.ipynb

例：

SMILES

→ RDKit 3D conformer

→ ASE Atoms

→ optional MACE optimization

→ SDF / XYZ保存


## 02_parameterization_review.ipynb

RESP/GAFF parameterを読み込み、

- atom
- atom type
- charge
- total charge
- parameter summary

を確認する。

parameterizationそのものをNotebook依存にはしない。


## 03_trajectory_analysis.ipynb

MDAnalysis等を使用し、

- trajectory load
- RMSD等の基本確認
- RDF
- energy data読み込み

の例を提供する。


# 25. CLI / workflow

v1では複雑なworkflow engineを導入しない。

Snakemake等は将来拡張とする。

基本操作例：

```bash
micromamba activate gmx-md

python scripts/prepare/build_system.py \
    --config systems/example/system.yaml

python scripts/prepare/parameterize.py \
    --config systems/example/system.yaml

./scripts/run/em.sh SYSTEM=example
./scripts/run/nvt.sh SYSTEM=example
./scripts/run/npt.sh SYSTEM=example
./scripts/run/production.sh SYSTEM=example
```

実際のCLI形式は実装時に多少改善してよい。


# 26. Makefile

以下のshortcutを用意することが望ましい。

```bash
make doctor
make test
make notebook
make probe
```

ただしMakefileへ重要なロジックを直接書かない。

処理本体はPythonまたはShellへ置く。


# 27. README

READMEには最低限以下を書く。

1. Repository目的
2. Architecture
3. Mac初期セットアップ
4. Micromambaセットアップ
5. environment作成
6. MACEセットアップ
7. GROMACS確認
8. Jupyter起動
9. example workflow
10. Linux server移行手順
11. GAFF/RESP workflow
12. Git管理方針
13. troubleshooting


# 28. MacからLinux serverへの移行

基本workflow：

```text
Mac
 ↓
git commit
 ↓
git push
 ↓
Linux server
 ↓
git clone / pull
 ↓
probe_server.sh
 ↓
AIによるserver環境判断
 ↓
Micromamba environment構築
 ↓
GROMACS設定
 ↓
SLURM設定
 ↓
doctor.sh
 ↓
small test
 ↓
production MD
```


# 29. 安全性

Setup scriptは以下を禁止する。

- sudoの無断実行
- system Pythonの変更
- system packageの削除
- shell設定ファイルの破壊的変更
- CUDA driverの自動変更
- NVIDIA driverの自動変更
- existing GROMACSの上書き
- cluster moduleの変更
- user data削除

破壊的操作が必要な場合は実行せず、必要なcommandをユーザーへ提示する。


# 30. v1の対象外

以下はv1では実装しなくてよい。

- 複雑な高分子自動構築
- polymer repeat-unit自動parameterization
- 無機結晶
- MACE-GROMACS直接連携
- QM/MM
- enhanced sampling
- free-energy calculation
- Snakemake
- Nextflow
- Kubernetes
- production Docker container
- automatic cloud deployment

ただし将来追加できる構造にしておく。


# 31. 完成条件

以下をすべて満たせばv1完成とする。

### Environment

- Micromamba環境を作成できる
- Jupyter kernelとして認識できる
- RDKitがimportできる
- ASEがimportできる
- MDAnalysisがimportできる
- AmberToolsが利用できる
- ACPYPEが利用できる
- MACEがimportできる


### GROMACS

- 外部GROMACSを検出できる
- `gmx --version` を記録できる
- GROMACSがなくてもPython環境は利用できる


### Workflow

- RDKitで3D分子を作成できる
- ASEへ渡せる
- optional MACE計算ができる
- RESP/GAFF parameterization pipelineの入口が存在する
- AMBER→GROMACS conversion pipelineが存在する
- MDP templateが存在する
- EM/NVT/NPT/production scriptが存在する


### Reproducibility

- Git commitを記録できる
- environment情報を記録できる
- input hashを記録できる
- parameterization条件を記録できる


### Repository

- trajectoryがGitへ入らない
-仮想環境がGitへ入らない
- scientific inputはGit管理可能
- READMEだけで新しいマシンへの導入手順が理解できる


# 32. 設計思想

本repositoryの最重要原則は以下。

```text
Notebook = 考える場所
Python = 再現する場所
Shell/SLURM = 計算する場所
Git = 条件を残す場所
GROMACS = classical MD engine
MACE = 補助的ML potential tool
Micromamba = Python環境管理
```

MacとLinuxを完全に同一環境にすることを目的としない。

代わりに、

「科学的条件は同じ、実行環境はマシンに最適化」

を基本思想とする。