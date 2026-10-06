# Ethanol 27-Molecule NVT Benchmark Tutorial 仕様書 v1.0

## 1. 目的

本チュートリアルでは、ASE + MACEを用いて実行した純エタノール27分子の短時間NVT-MD計算に対し、GROMACS + 古典力場を用いて、**力場以外の条件を可能な限り対応させた比較計算**を構築する。

本計算の目的は、

- MACE-MDと古典力場MDの完全一致
- 平衡液体エタノールの高精度物性評価

ではない。

目的は、

> 同一に近い初期構造・初期速度・温度制御・時間刻み・計算時間・解析条件を与えたとき、MACE potentialと古典力場でどのような違いが生じるかを理解する教材計算を作ること

である。

特に、

```text
ASE / MACE
        ↓
同じ構造・MD条件
        ↓
GROMACS / GAFF + RESP
        ↓
同一のPython RDF解析
```

という比較可能なworkflowを構築する。


# 2. 前提環境

本チュートリアルは、別途構築済みの

`GROMACS MD Workbench`

上で実行する。

想定環境：

- ASE 3.29.0
- mace-torch 0.3.16
- NumPy
- RDKit
- MDAnalysis
- AmberTools
- ACPYPE
- GROMACS
- Jupyter
- Python

Python環境はMicromambaで管理する。

GROMACSは外部commandとして利用する。

本チュートリアルのために、ASE 3.29.0およびmace-torch 0.3.16を再現可能な環境情報を残す。


# 3. 比較の基本原則

## 3.1 対応させる条件

可能な限り以下をMACE版とGROMACS版で一致させる。

- 化学組成
- 分子数
- 原子数
- 初期分子構造
- 初期配置
- periodic box
- 初期密度
- 初期速度
- 乱数seed
- COM momentum除去
- 初期温度
- constraintなし
- timestep
- thermostatの考え方
- thermostat chain length
- equilibration時間
- production時間
- trajectory保存間隔
- RDF解析frame
- RDF binning
- RDF finite-size normalization


## 3.2 一致させられない条件

potential energy surfaceは異なる。

### 元計算

```text
MACE-OFF
mace_off(model="small",
         device="cpu",
         default_dtype="float64")
```

### GROMACS計算

```text
GAFF / GAFF2
+
RESP charge
+
classical bonded/nonbonded interaction
```

したがって、trajectory、energy、temperature fluctuation、RDF等が数値的に一致することを要求しない。


# 4. 系の定義

## 4.1 分子

純エタノール：

```text
CH3CH2OH
```

分子数：

```text
27 molecules
```

1分子：

```text
9 atoms
```

総原子数：

```text
27 × 9 = 243 atoms
```

追加の、

- solvent
- ion
- counter ion

は使用しない。


# 5. 初期エタノール構造

ASEを使用して、

```python
ase.build.molecule("CH3CH2OH")
```

からエタノール1分子を生成する。

その分子について、

1. center of massを原点へ移動
2. 分子のx軸を `(1, 1, 1)` 方向へ向ける

というorientation処理を行う。

この処理を再現可能なPython関数として実装する。


# 6. 27分子の配置

## 6.1 格子

```text
3 × 3 × 3
```

の格子点へエタノールを配置する。

分子間隔：

```text
4.6 Å
```

すなわち、

```text
0.46 nm
```

とする。


## 6.2 ランダム微小回転

各分子へ、

```text
x axis
y axis
z axis
```

それぞれについて、

```text
-5° ～ +5°
```

の一様乱数回転を加える。

乱数生成器：

```python
numpy.random.default_rng(42)
```

を使用する。

Python global RNG等へ変更しない。


# 7. Simulation box

立方体セル：

```text
13.8 × 13.8 × 13.8 Å
```

すなわち、

```text
1.38 × 1.38 × 1.38 nm
```

とする。

periodic boundary condition：

```text
x = True
y = True
z = True
```

box volumeはMD中固定する。

barostatは使用しない。


# 8. 初期密度

初期配置の質量密度：

```text
0.786 g/cm³
```

を基準値とする。

構造生成後、実際の、

- atom count
- molecular mass
- box volume
- calculated density

をPythonで計算し、reportへ出力する。

期待値と大きく異なる場合はwarningを出す。


# 9. Energy Minimization

元MACE計算では、初期構造に対するenergy minimizationを実施していない。

したがってGROMACS版でも、

```text
Energy Minimizationなし
```

を基本条件とする。

ただし、古典力場で初期構造を評価した際に、

- 異常な原子重なり
- 非有限energy
- 明らかな巨大force
- topology error
- 不正なbond
- catastrophic clash

などが存在しないことをpreflight checkする。

問題が検出された場合、Codexまたはscriptが勝手にenergy minimizationを実行して条件を変更してはならない。

その場合は計算を停止し、

```text
指定された初期条件ではGROMACS古典力場計算を
安全に開始できない可能性がある
```

とreportする。


# 10. Classical force field

## 10.1 方針

GROMACS版では、

```text
GAFF または GAFF2
+
RESP charge
```

を使用することを想定する。

使用したforce fieldを明示的に記録する。


## 10.2 必須記録項目

最低限、以下をdocumentationへ残す。

- GAFF / GAFF2
- atom type
- atomic charge
- total molecular charge
- bond parameter
- angle parameter
- dihedral parameter
- improper parameter
- Lennard-Jones parameter
- combination rule
- exclusion rule
- 1–4 Coulomb scaling
- 1–4 LJ scaling
- electrostatics method
- Coulomb cutoff
- van der Waals cutoff
- long-range treatment


# 11. RESP charge

RESP chargeを使用する。

以下についてはscientific inputとし、自動推測で変更しない。

- QM backend
- QM method
- basis set
- molecular geometry
- total charge
- multiplicity
- ESP calculation condition
- RESP fitting condition

これらはconfiguration fileへ記録する。

RESP情報が存在しない場合、

```text
AM1-BCCへ自動fallback
```

してはならない。

必要情報不足として停止する。


# 12. Constraints

元計算ではconstraintを使用していない。

したがってGROMACS版でも、

```text
constraints = none
```

とする。

特に、

- all-bonds
- h-bonds
- rigid ethanol

は使用しない。

全243原子を明示的に積分する。


# 13. Initial velocities

## 13.1 温度

```text
300 K
```

のMaxwell–Boltzmann分布から初期運動量を生成する。


## 13.2 RNG

元計算と同様、

```python
numpy.random.default_rng(42)
```

を使用する。


## 13.3 COM removal

初期速度生成後に、

```text
system全体のcenter-of-mass translational momentum
```

を1回だけ除去する。

各エタノール分子ごとの、

- translational momentum
- rotational momentum

は除去しない。


## 13.4 Temperature rescaling

COM removal後、

瞬間温度が

```text
300 K
```

ちょうどになるよう全原子のmomentumを一様に1回だけscaleする。

このscaleはinitialization時のみである。

MD中にvelocity rescaling thermostatを使用しない。


## 13.5 GROMACSへの受け渡し

GROMACS自身に別のrandom velocityを生成させない。

したがって、

```ini
gen_vel = no
```

とする。

Python側で生成した速度をGROMACS inputへ保持して渡す。

これにより、

単にrandom seedを一致させるのではなく、

```text
initial positions
+
initial velocities
```

そのものを可能な限り一致させる。


# 14. Degrees of freedom

全原子を可動とする。

constraintなし。

COM removalによる自由度補正を含め、

初期温度計算に使用した自由度の定義をコードまたはmetadataに記録する。

ASE側とGROMACS側で自由度の定義が異なる場合は明記する。


# 15. MD ensemble

```text
NVT
```

とする。

- particle number: fixed
- volume: fixed
- temperature: 300 K

pressure controlは使用しない。


# 16. Integrator

元計算のNosé–Hoover chain条件へ対応させるため、GROMACSではvelocity-Verlet系integratorを使用する。

基本候補：

```ini
integrator = md-vv
```

とする。

leap-frog `md` へ無条件に変更しない。


# 17. Timestep

```text
0.5 fs
```

すなわち、

```text
0.0005 ps
```

とする。

GROMACS：

```ini
dt = 0.0005
```

とする。


# 18. Thermostat

## 18.1 元条件

ASE：

```text
Nosé–Hoover chain
chain length = 3
tdamp = 50 fs
```

## 18.2 GROMACS

thermostat：

```ini
tcoupl = nose-hoover
```

chain length：

```ini
nh-chain-length = 3
```

を使用する。


## 18.3 Damping timeの対応

ASEの `tdamp = 50 fs` をGROMACS `tau_t = 0.05 ps` と単純対応させない。

thermostat massの定義差を考慮する。

thermostat massを対応させる場合、

```text
tau_t = 2π × tdamp
```

とし、

```text
tdamp = 0.05 ps
```

から、

```text
tau_t ≈ 0.314159 ps
```

を候補値とする。

最終的なmappingの根拠をdocumentationへ記載する。


# 19. Temperature coupling group

全エタノール原子を1つのtemperature coupling groupとして扱う。

分子ごとのthermostat couplingは行わない。


# 20. Pressure coupling

なし。

```ini
pcoupl = no
```

とする。


# 21. Simulation duration

## 21.1 Equilibration

```text
200 steps
```

timestep：

```text
0.5 fs
```

なので、

```text
100 fs
=
0.1 ps
```

とする。


## 21.2 Production

続く、

```text
400 steps
```

すなわち、

```text
200 fs
=
0.2 ps
```

をproduction区間とする。


## 21.3 Total

合計：

```text
600 steps
300 fs
0.3 ps
```


# 22. 1本のcontinuous simulationとして実行

EquilibrationとProductionを別々のMD runへ分割しない。

原則、

```text
step 0
↓
step 600
```

まで1本の連続NVT計算として実行する。

理由：

- thermostat内部状態の連続性を維持
- restart差を避ける
- 元ASE計算との対応を簡潔にする

解析時のみ、

```text
step <= 200
```

をequilibration、

```text
step > 200
```

をproductionとして扱う。


# 23. Trajectory保存

保存間隔：

```text
5 steps
```

とする。

時間では、

```text
5 × 0.5 fs
=
2.5 fs
=
0.0025 ps
```

である。

圧縮trajectoryを使用してよい。

例：

```ini
nstxout-compressed = 5
```

ただしRDF精度へ影響する座標quantization設定を確認する。


# 24. Production frames

解析対象は、

```text
step > 200
```

のみ。

対象step：

```text
205
210
215
...
595
600
```

総frame数：

```text
80 frames
```

対象時刻：

```text
102.5 fs ～ 300.0 fs
```

すなわち、

```text
0.1025 ps ～ 0.3000 ps
```


# 25. RDF解析

標準GROMACS RDFだけに依存せず、

Python + NumPy / MDAnalysisで、

元MACE trajectoryとGROMACS trajectoryの両方へ適用可能なcustom RDF実装を作る。


# 26. RDF種類

以下の2種類を計算する。

## 26.1 O–O RDF

各エタノール分子のhydroxyl oxygen同士。

```text
O_i → O_j
```

ただし、

```text
i ≠ j
```

とする。


## 26.2 O–H RDF

中心：

```text
hydroxyl oxygen
```

相手：

```text
別分子のhydroxyl hydrogen
```

とする。

同一エタノール分子中のO–H pairは除外する。

つまり、

```text
O_i → H_j
i ≠ j
```

のみを数える。


# 27. Periodic distance

周期境界条件を適用し、

```text
minimum image convention
```

による最短距離を使用する。


# 28. RDF範囲

box length：

```text
13.8 Å
```

に対し、

```text
r_max = 0.49 × 13.8 Å
```

とする。

したがって、

```text
r_max = 6.762 Å
```

である。


# 29. RDF bins

```text
45 bins
```

とする。

bin width：

```text
6.762 / 45
≈ 0.15027 Å
```

すなわち約、

```text
0.1503 Å
```

である。


# 30. RDF normalization

各中心原子について、各binのpair countを、

```text
spherical shell volume
×
partner number density
```

で規格化する。


## 30.1 partner number density

27分子から自己分子を除外するので、

```text
Npartner = 26
```

とする。

したがって、

```text
ρpartner = 26 / V
```

を使用する。

ここで `V` はsimulation box volume。


## 30.2 finite-size correction

通常のbulk RDFで使われる、

```text
27 / V
```

ではなく、

```text
26 / V
```

を明示的に使用する。

これは小さい閉鎖系に対する有限個数補正として扱う。

GROMACSやMDAnalysisの標準RDF normalizationと結果が一致することを要求しない。


# 31. RDF実装の確認

custom RDF implementationについて、最低限以下をunit testする。

- same-molecule pairが除外されている
- periodic boundaryが適用される
- r_maxを超えるpairを数えない
- 45 binsになる
- 26/V normalizationになる
- frame数が80になる
- step 200を含まない
- step 205を含む
- step 600を含む


# 32. MACE版trajectoryとの共通解析

可能であれば、RDF解析ロジックは、

```text
MACE trajectory
GROMACS trajectory
```

両方に同一関数を使用する。

trajectory readerのみadapter化してよい。

解析本体では、

```text
coordinates
box
molecule IDs
O atom indices
hydroxyl-H indices
```

を共通形式へ変換する。


# 33. Directory構造

本チュートリアルは例えば以下へ配置する。

```text
systems/
└── ethanol_27/
    ├── README.md
    ├── system.yaml
    │
    ├── structures/
    │   ├── ethanol_single.*
    │   ├── ethanol27_initial.*
    │   └── ethanol27_with_velocities.*
    │
    ├── parameters/
    │   ├── resp/
    │   ├── amber/
    │   └── gromacs/
    │
    └── analysis/
        └── reference/

examples/
└── ethanol_nvt_comparison/
    ├── README.md
    ├── build_initial_structure.py
    ├── generate_velocities.py
    ├── preflight.py
    ├── run_nvt.sh
    ├── analyze_rdf.py
    └── notebooks/
        └── ethanol_comparison.ipynb
```

実際のrepository architectureへ合わせ、過剰な重複は避ける。


# 34. system.yaml

最低限以下を記述できるようにする。

```yaml
system:
  name: ethanol_27
  molecule: ethanol
  molecule_count: 27
  atoms_per_molecule: 9

box:
  type: cubic
  length_A: 13.8
  pbc: [true, true, true]

placement:
  grid: [3, 3, 3]
  spacing_A: 4.6
  rotation_deg:
    min: -5.0
    max: 5.0
  rng_seed: 42

velocities:
  temperature_K: 300.0
  rng_seed: 42
  remove_system_com: true
  remove_molecular_com: false
  remove_molecular_rotation: false
  exact_temperature_rescale: true

force_field:
  family: gaff2
  charge_model: resp

md:
  ensemble: NVT
  temperature_K: 300.0
  timestep_fs: 0.5
  steps: 600

thermostat:
  type: nose-hoover-chain
  chain_length: 3
  ase_tdamp_fs: 50.0
  gromacs_tau_t_ps: 0.3141592654

trajectory:
  stride_steps: 5

analysis:
  production_start_exclusive_step: 200
  r_max_A: 6.762
  rdf_bins: 45
  finite_size_normalization: true
```


# 35. GROMACS MDP

MDPは本チュートリアル専用fileとして管理する。

最低限明示：

```text
integrator
dt
nsteps
constraints
tcoupl
tau_t
ref_t
nh-chain-length
pcoupl
nstxout-compressed
nstenergy
nstlog
periodic boundary
electrostatics
vdW
cutoff
```

force fieldに関係する、

- cutoff
- PME等のelectrostatics method
- dispersion correction

については、使用したGAFF/RESPモデルとの整合性を確認して設定する。

元MACE計算に存在しない条件なので、documentationへ明記する。


# 36. Metadata

各runについて最低限以下を保存する。

```text
ASE version
mace-torch version
NumPy version
GROMACS version
AmberTools version
ACPYPE version
GAFF/GAFF2
RESP settings
initial structure hash
initial velocity hash
MDP hash
topology hash
git commit
random seeds
box
density
timestep
nsteps
thermostat
tau_t
chain length
constraints
electrostatics method
cutoffs
start time
command
```


# 37. Validation

simulation開始前に以下を確認する。

## Structure

- 27 molecules
- 243 atoms
- 9 atoms/molecule
- box = 1.38 nm cubic
- PBC = xyz
- density ≈ 0.786 g/cm³


## Velocity

- RNG seed = 42
- system COM ≈ 0
- instantaneous temperature ≈ 300 K
- molecular COMは除去されていない


## Topology

- 27 ethanol molecules
- total charge整合
- constraintsなし
- missing parameterなし


## MD

- dt = 0.0005 ps
- nsteps = 600
- NVT
- target = 300 K
- Nose–Hoover
- NHC length = 3
- no barostat
- no generated GROMACS velocity


# 38. Output

最低限以下を生成する。

```text
initial structure
initial velocities
topology
MDP
TPR
trajectory
energy
final structure
metadata
temperature trace
potential energy trace
RDF O–O
RDF O–H
```


# 39. Visualization

Notebookでは最低限以下を表示する。

### Initial structure

- simulation box
- 27 ethanol molecules

### MD

- Temperature vs time
- Potential Energy vs time

### RDF

同一グラフ上に必要に応じて、

```text
MACE O–O
GROMACS O–O
```

および、

```text
MACE O–H
GROMACS O–H
```

を比較表示する。

比較時は、

```text
Potentialが異なるため一致を期待していない
```

ことを明記する。


# 40. Interpretation

本simulationは、

```text
equilibration = 100 fs
production = 200 fs
```

しかない。

したがって、

- equilibrium liquid structure
- converged RDF
- experimental RDF
- diffusion coefficient
- reliable thermodynamic property

を得たとは主張しない。

元MACE計算同様、初期3×3×3格子配置の影響がproductionにも残る可能性が高い。

結果は、

> 短時間MDにおけるMACEとGAFF/RESPの挙動差を理解する教材

として解釈する。


# 41. 本チュートリアルで答えたい問い

最終的なNotebook / READMEでは、少なくとも以下を考察できるようにする。

1. 同じ初期座標・初期速度から開始しても、MACEとGAFF/RESPでtrajectoryがなぜ急速に異なるのか。
2. thermostat条件を対応させても完全に同じNVT dynamicsにならない可能性があるのはなぜか。
3. constraintなしで0.5 fs timestepを使う理由は何か。
4. 短いproductionではRDFが収束しないのはなぜか。
5. O–O RDFとO–H RDFはそれぞれ何を表すか。
6. 初期格子構造がRDFへどのように現れるか。
7. 26/V normalizationと27/V normalizationで何が違うか。
8. 古典力場とML potentialの比較で、「条件を揃える」とは何を意味するか。


# 42. 実装時にAIが勝手に変更してはいけない条件

以下は明示的な指示なしに変更しない。

```text
27 molecules
243 atoms
13.8 Å box
3×3×3 placement
4.6 Å spacing
rotation seed = 42
velocity seed = 42
300 K
no minimization
no constraints
0.5 fs timestep
NVT
NHC chain length = 3
ASE tdamp = 50 fs
600 steps
trajectory every 5 steps
production = step > 200
r_max = 6.762 Å
45 bins
26/V normalization
```


# 43. AIが判断してよい項目

以下については、既存Workbench architectureやsoftware制約を見て適切な方法を選択してよい。

- file format
- internal function names
- Python class設計
- command line interface
- directory細部
- logging implementation
- unit test framework細部
- plotting implementation
- trajectory adapter
- temporary file management

ただし科学的条件を変えない。


# 44. 未確定の科学的項目

以下は本仕様書時点で最終決定していない。

- GAFF or GAFF2
- RESP生成に使用するQM backend
- RESP QM method / basis
- exact electrostatic treatment
- cutoff values
- dispersion correction
- exact 1–4 settings

これらはforce field整合性を確認して決定し、採用理由をREADMEへ記録する。

AIがsilent defaultで決めてはならない。


# 45. 完成条件

本チュートリアルは以下を満たした場合に完成とする。

### Build

- ASEから27分子構造を再現できる
- random rotationsがseed 42で再現できる
- box/densityが仕様通りである

### Velocity

- NumPy RNGで初期速度を生成できる
- system COM momentumを除去できる
- 300 Kへ1回だけscaleできる
- GROMACSへ同じvelocityを渡せる

### Force field

- RESP/GAFF系topologyが作成されている
- parameter情報が記録されている
- constraintなしである

### MD

- 600-step continuous NVTが実行できる
- NHC length 3を使用できる
- 0.5 fs timestepである
- 5-stepごとにtrajectoryを保存できる

### Analysis

- production frameが80 frameである
- O–O RDFを計算できる
- intermolecular O–H RDFを計算できる
- 26/V normalizationを実装している
- MACE/GROMACS双方へ同じRDF logicを適用できる

### Documentation

- MACEと古典力場の差を明記している
- thermostat mappingを説明している
- 短時間simulationの限界を明記している
- 未確定条件をsilent defaultにしていない


# 46. 最終的な教材workflow

```text
ASE molecule("CH3CH2OH")
          ↓
single ethanol orientation
          ↓
3 × 3 × 3 placement
          ↓
random rotation
seed = 42
          ↓
13.8 Å cubic box
          ↓
RESP + GAFF/GAFF2
          ↓
NumPy initial velocities
seed = 42
          ↓
COM removal
          ↓
exactly 300 K scaling
          ↓
GROMACS preflight
          ↓
600-step NVT
dt = 0.5 fs
          ↓
trajectory every 5 steps
          ↓
step > 200
          ↓
80 production frames
          ↓
custom RDF
O–O / O–H
26/V normalization
          ↓
MACE resultとの比較
          ↓
教育的考察
```


# 47. 設計思想

本チュートリアルの最重要原則は、

> Potential以外の条件を可能な限り揃え、揃えられない条件は明示する。

ことである。

「GROMACSでMACEを再現する」のではなく、

> MACEによるMDと古典力場によるMDを、可能な限り公平な条件で比較する

ことを目的とする。

また、本計算は平衡液体エタノールの本格simulationではなく、

```text
MD workflow
thermostat
force field
initial conditions
RDF
finite-size effect
ML potential vs classical potential
```

を学習するための短時間教材計算として位置付ける。