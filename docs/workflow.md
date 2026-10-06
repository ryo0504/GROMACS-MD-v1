# Workflow と MDP

SMILES/SDF/MOL2 → RDKit explicit-H構造 → ASE → optional MACE → external QM ESP/RESP → GAFF/GAFF2 → AMBER → ACPYPE → periodic system → EM → NVT → NPT → production → MDAnalysis。
MACE利用時はmodel file/device/fmax/stepsをconfigへ明示します。modelのhashを記録し、自動model選択/downloadはしません。
MACE最適化後にESPを計算し、RESPと一致したgeometry/atom順序を保持してください。

## 教育用MDPの意味

MDPはperiodic condensed-phaseを想定した出発点です。気相単分子や希薄系、全有機物に共通の推奨値ではありません。
unitsは時間ps、距離nm、temperature K、pressure barです。GRO nmに対しASE/RDKit/MDAnalysisの座標は通常Åです。

| parameter | 意味とこのテンプレート |
|---|---|
| integrator | steepはEM、mdはleap-frog MD |
| nsteps | 最大step数。EM5000、NVT/NPT50000、MD500000 |
| dt | MD timestep 0.002ps。制約/振動数/力場に合わせて検証 |
| emtol, emstep | EMのforce停止基準と最大step長 |
| continuation | checkpoint/前段geometryによる継続処理 |
| constraints | MDではh-bonds制約。EMはnone |
| constraint-algorithm, lincs-order | LINCS constraint solverと展開次数 |
| cutoff-scheme, nstlist, rlist | Verlet neighbor list、更新周期、list距離 |
| coulombtype, rcoulomb | PME electrostaticsとreal-space cutoff |
| vdwtype, rvdw | LJ cutoff。1.0nmが適切か対象力場で確認 |
| DispCorr | uniform-fluid仮定のLJ長距離energy/pressure補正 |
| pbc | xyz方向周期境界 |
| tcoupl, tc-grps | v-rescale thermostat、System全体でcoupling |
| tau-t, ref-t | 温度coupling時間、温度。ref-tはconfigから反映 |
| pcoupl, pcoupltype | NVTなし、NPT C-rescale、MD Parrinello-Rahman、isotropic |
| tau-p, ref-p | pressure relaxation timeと目標pressure。ref-pはconfigから反映 |
| compressibility | 4.5e-5/barは水に近い例。対象有機物に対して利用者が決め直す |
| gen-vel, gen-temp, gen-seed | NVTのみ初期velocity生成。temperature/seedはconfigから反映 |
| nstenergy, nstlog, nstxout-compressed | energy/log/compressed trajectory出力周期 |

algorithm詳細は [GROMACS MDP reference](https://manual.gromacs.org/current/user-guide/mdp-options.html) を確認してください。
C-rescaleはGROMACS 2021以降です。古いserver版では条件を勝手にfallbackさせずMDPを明示的に見直してください。

## Runの記録

runごとに新しいdirectoryを用意します。`prepare-run` が科学入力とinclude graphを保存し、後からの入力改変をhashで検出します。
gromppは各stageにprocessed.mdp/processed.topを出力し、warningを無視しません。
NVTはEM結果から開始しvelocityを生成、NPTはNVT checkpoint、productionはNPT checkpointを受け取ります。
途中中断の自動resumeはv1では行いません。新しいrunで再実行するか、利用者がcheckpointとappend方針を明示してGROMACSで復旧し、metadataを保存してください。
stage失敗時はmetadata command statusとgrompp.log/mdrun.logを確認してください。

## 解析

Notebook 03はRMSD、RDFと `gmx energy` でexportしたXVGを読みます。
RMSD前にはPBC補正/unwrap/centeringを行い、reference/selectionとfitの意味を確認してください。
RDFはÅ単位で、intramolecular/self pair exclusionsとbulk density normalizationの科学的判断が必要です。
