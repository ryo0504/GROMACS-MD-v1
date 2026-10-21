# GROMACS MD Workbench

Macで入力作成と小規模確認を行い、Linux/HPCで本番の classical MD を実行する作業環境です。
Notebook は探索・解析、Python は再現可能な処理、Shell/SLURM は計算実行を担います。
Micromamba `gmx-md` に科学 Python をまとめ、GROMACS はホスト側の外部 command とします。
MACE は構造最適化・energy/force評価用の補助ツールです。

## Macで開始

このPCの構築結果と実際の環境パスは [環境構築記録](docs/current_pc_setup.md) を参照してください。Notebookのカーネルは **Python (GROMACS MD Workbench)**（ID: `gmx-md`）です。

このディレクトリで次を実行すると、構築済み環境でJupyterLabを開けます。

```bash
./scripts/setup/start_jupyter.sh
```

terminalでもPython/GROMACSを使う場合：

```bash
source scripts/setup/common.sh
eval "$("$MM" shell hook -s zsh)"
micromamba activate gmx-md
python --version
gmx --version
make doctor
```

このPCではbinary environmentを同期対象外の `~/.local/share/mdworkbench/GROMACS-MD-v1/` に配置しています。`.local/host.env` からこのPCの設定を読み込みます。`MAMBA_ROOT_PREFIX` や `GMX_BIN` をterminalで明示指定した場合は、その設定を優先します。

別のPCで新規構築するときは `scripts/setup/setup_micromamba.sh`、`scripts/setup/install_mace.sh`、`scripts/setup/install_gromacs_mac.sh` を使い、そのマシンのOS/architectureに対応させてください。Mac/Linuxでbinary environmentをそのままコピーせず、環境定義とlockから再構成します。

## Example workflow

`systems/example/system.yaml` は **ethanolの教育用条件テンプレート** です。RESP source は未指定なのでparameterizationは意図どおり停止します。
`defaults.yaml` は参照専用で、自動的に科学条件へmergeしません。相対パスは YAML の位置から解決します。

```bash
python scripts/prepare/build_system.py \
  --config systems/example/system.yaml --output runs/example-build
```

同じseedで生成した教育用 SDF/XYZ とconstruction情報を `systems/example/structures/` に保存済みです。
実際のQM/RESPには適切なgeometryを確認・保存し、外部 QM/RESP 結果を用意します。
`system.yaml` の QM method/basis、charge、multiplicity、source、atom order、fitting条件を実際の計算に合わせて編集してください。

```bash
python scripts/prepare/parameterize.py \
  --config systems/example/system.yaml --output runs/example-parameters
```

この処理は RESP fit/import → antechamber (`-c rc`) → parmchk2 → tleap → ACPYPE (`-p`, `-x`) を実行します。
AM1-BCCへのfallbackはありません。最終 MOL2、frcmod、RESP input/charge、metadataとtop/itpなどの小さい確定入力を `systems/` に保存してGit管理できます。

次に物理的に妥当なbox、分子数、組成、溶媒、densityをもつsystem.gro/system.topを利用者が構築します。
単分子のACPYPE出力はbulk系ではありません。v1のbuild_systemは単分子生成用で、自動packing/solvationは含みません。

```bash
mdworkbench prepare-run --config systems/example/system.yaml \
  --structure systems/example/structures/system.gro \
  --topology systems/example/parameters/system.top \
  --mdp-dir mdp --output runs/example-001
./scripts/run/em.sh --run-dir runs/example-001
./scripts/run/nvt.sh --run-dir runs/example-001
./scripts/run/npt.sh --run-dir runs/example-001
./scripts/run/production.sh --run-dir runs/example-001
```

`prepare-run` はatom/residue数・atom順序・total charge・box・includeを確認して入力をsnapshotします。
includeが別の力場ディレクトリにある場合は `--include-dir /path/to/gromacs/share/gromacs/top` を追加します。
temperature/pressure/seedはsystem.yamlからコピーしたMDPへ明示的に反映し、dt/cutoff/thermostat/barostat等は選択したMDPに従います。
gromppのwarningは無視しません。専用600-stepエタノール教材の既知COM警告だけは、[内容照合付きの限定例外](examples/ethanol_nvt_comparison/README.md#入力が揃った後の操作)があります。NPT/productionへは前段のcheckpointを渡します。
同時実行と出力上書きは拒否します。失敗runのログを確認し、修正した入力で新しいrunを準備してください。

## Linux server

```bash
git clone <YOUR_REPOSITORY_URL>
cd <CLONED_REPOSITORY>
./scripts/setup/probe_server.sh
cat cluster/detected/platform_report.txt
```

reportとcluster policyから既存module/GROMACS、CUDA、MPI、SLURM資源を選び、[server_setup](docs/server_setup.md)の手順で進めます。
server固有値は `config/platform.local.yaml` に保存し、run scriptへ `--platform` で明示的に渡します。
SLURMテンプレートのplaceholderを埋める前にsubmitしないでください。

## 再現性とGit

各処理は `metadata.json` にinput SHA256、config、git commit/dirty、host、packages、GROMACS build、GPU/CUDAなどを記録します。
runはexecuted commandとstatus、gromppのprocessed topology/MDP/TPR hashも記録します。
Condaのplatform別explicit lockとpip inventoryは `env/locks/` へexportできます。

```bash
./scripts/setup/export_environment.sh
```

`runs/`、trajectory/checkpoint、cache、環境、ML modelはignoredです。小さいMDP/TOP/ITP/MOL2/GRO/RESP条件は管理可能です。
Gitはmain branchで初期化済みです。最初のcommitとremote設定は利用者のリモート運用方針に合わせて行ってください。
未commit状態のmetadataではcommitがnullになり、偽のcommitを記録しません。

詳しくは [architecture](docs/architecture.md)、[workflow/MDP](docs/workflow.md)、[parameterization](docs/parameterization.md)、[environment](env/README.md)、[troubleshooting](docs/troubleshooting.md) を参照してください。
実機での結果と検証の範囲は [v1検証記録](docs/validation.md) に保存しました。

## 確認した公式資料

- [Micromamba installation](https://mamba.readthedocs.io/en/stable/installation/micromamba-installation.html)
- [AmberTools](https://ambermd.org/AmberTools.php) / [Amber manual](https://ambermd.org/doc12/Amber25.pdf)
- [Antechamber official source/CLI](https://github.com/Amber-MD/AmberClassic/tree/main/src/antechamber)
- [ACPYPE conversion](https://acpype.readthedocs.io/en/latest/)
- [MACE installation](https://mace-docs.readthedocs.io/en/latest/guide/installation.html)
- [GROMACS installation](https://manual.gromacs.org/current/install-guide/index.html) / [MDP reference](https://manual.gromacs.org/current/user-guide/mdp-options.html)
- [Homebrew GROMACS](https://formulae.brew.sh/formula/gromacs)

## エタノール27分子のGROMACS習熟教材

仕様書に対応する [Ethanol NVT comparison tutorial](examples/ethanol_nvt_comparison/README.md) と [Notebook](examples/ethanol_nvt_comparison/notebooks/ethanol_comparison.ipynb) があります。構造生成→速度→RDF演習→RESP/GAFF点検→preflight→600-step連続NVT→結果の読み方を学べます。第5節には無料のPySCFによるQM/ESP→2段階RESP→GAFF2の実行手順を追加しました。RESP等の科学入力が未準備でも前半の演習は実行できます。汎用EM/NVT/NPT workflowとは別の専用教材です。
