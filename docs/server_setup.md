# Linux server adaptation

clone後、最初に `./scripts/setup/probe_server.sh` を実行し `cluster/detected/platform_report.txt` を読みます。
login shellでのみmoduleが初期化される環境では `module list` / `module avail` も確認しreportへ追記してください。
reportはCPU/RAM、GPU/driver、nvcc、compiler/CMake、MPI/OpenMPI、SLURM、filesystem/quotaを可能な範囲で報告します。
module loadやdriver変更、sudoは自動実行しません。

1. cluster policyに従って既存GROMACS moduleを選ぶ。
2. `gmx --version` または `gmx_mpi --version` からGPU/MPI/buildを確認。
3. single-node thread-MPIかexternal MPIかを選択。external MPIならsite推奨のsrun/mpirunをlauncherに設定。
4. Micromamba環境を作成し、hostに適切なMACE/PyTorch backendを選ぶ。
5. SLURM templateのpartition/account/CPU/GPU/memory/timeを埋め、log directoryを事前に作成。
6. doctor/pytestとsmall testで確認後、productionをsubmit。

```bash
./scripts/setup/setup_micromamba.sh
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
export PATH="$PWD/.local/bin:$PATH"
eval "$(micromamba shell hook -s bash)"
micromamba activate gmx-md
# module load <site-approved-gromacs> / source <GMXRC> as instructed by your site
./scripts/setup/detect_gromacs.sh
./scripts/setup/install_mace.sh
cp config/platform.example.yaml config/platform.local.yaml
# Edit host config and SLURM placeholders based on report/site policy.
make doctor
make test
```

Macのosx-arm64 lockをLinuxで使わず、Linux専用lockがなければYAMLからsolveしてexportしてください。
report中のfilesystem/hostname等は環境情報なので、外部に渡す場合は内容を確認してください。
