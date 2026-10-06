# エタノール27分子：GROMACS NVTを学ぶ

仕様書に対応する新しい教材です。まず **Notebookの第0〜4節を順番に実行**してください。初期構造・初期速度・RDFの演習はRESPやGROMACSがなくても実行できます。MDの結果は実際の出力があるときだけ表示します。

## 学習の道筋

| 段階 | 操作 | 自分で確認すること |
|---|---|---|
| 1 | ASEで1分子→27分子を作る | 配置・原子数・密度・周期境界 |
| 2 | Maxwell–Boltzmann運動量を生成 | COM除去と温度の一回のスケール |
| 3 | 初期配置のRDFを見る | O–H共有結合を数えない理由、26/Vの意味 |
| 4 | RESP/GAFFの入力を点検 | 電荷・原子順序・1–4相互作用 |
| 5 | `grompp`→ゼロステップ評価 | TPRとは何か、異常な力をどう検出するか |
| 6 | 600ステップを連続実行 | MDP、EDR、XTC、GRO、LOGの役割 |
| 7 | 時系列と80フレームのRDF | 平衡・収束を主張できない理由 |

教材Notebook：[ethanol_comparison.ipynb](notebooks/ethanol_comparison.ipynb)。処理本体は `src/mdworkbench/ethanol_benchmark.py` と `src/mdworkbench/analysis/ethanol_rdf.py` です。汎用WorkbenchのEM→NVT→NPT→Productionは使用しません。

## 最初の操作（repository rootで実行）

既存 `gmx-md` 環境が健全である場合：

```bash
source scripts/setup/common.sh
"$MM" run -n gmx-md python -m mdworkbench.ethanol_benchmark build
"$MM" run -n gmx-md jupyter lab examples/ethanol_nvt_comparison/notebooks/ethanol_comparison.ipynb
```

Notebookでは `Python (GROMACS MD Workbench)` 相当の `gmx-md` kernelを選びます。環境の作成・修復はrepositoryの `README.md` / `env/README.md` を参照してください。build先が既存なら新しい `--output runs/ethanol27_initial_02` を使います。既存runを上書きしません。Notebookの `BUILD` / `RUN` も対応させます。

このPCでは `gmx-md` を再構築し、**Python (GROMACS MD Workbench)** カーネルを登録しました。構築結果は [環境構築記録](../../docs/current_pc_setup.md) を参照してください。以下の `gmx-ethanol27` は別の学習用環境を作りたい場合の選択肢であり、このPCで追加作成する必要はありません（このenvironment solve自体は未検証）。

```bash
source scripts/setup/common.sh
"$MM" create -f env/ethanol27-study.yml
"$MM" run -n gmx-ethanol27 python -m pip install -e .
"$MM" run -n gmx-ethanol27 python -m ipykernel install --user --name gmx-ethanol27 --display-name "Python (Ethanol 27 tutorial)"
"$MM" run -n gmx-ethanol27 jupyter lab examples/ethanol_nvt_comparison/notebooks/ethanol_comparison.ipynb
```

この場合はNotebookのkernelを `Python (Ethanol 27 tutorial)` に変更し、以降のcommandも `-n gmx-ethanol27` を使います。MACEは比較データがないため必須dependencyに含めていません。将来MACEを実行する場合はPyTorch backendを確認してmace-torch0.3.16を別途明示的に導入し、実versionを記録します。GROMACSはこの環境に含めず外部commandを使います。

構造作成にはASE・RDKit・NumPy・PyYAML、Notebookにはmatplotlib、MD解析にはMDAnalysis・ParmEdが必要です。MACEの実行は不要です。buildのmetadataには実際のASE/MACE等のversionが保存されます。仕様のASE3.29.0/MACE0.3.16と異なる場合は再現の差として扱います。

## RESPとは何を用意することか

RESPはQM計算で得た分子周囲の静電ポテンシャルに原子電荷をフィットする方法です。SMILESだけからこの仕様のRESP電荷が決まるわけではありません。現在はQM/RESP入力が未準備です。MDを実行済みとは扱いません。

`systems/ethanol_27/system.yaml` を編集し、以下を埋めます。`null` は未決定であり、仮の電荷やAM1-BCCに置換されません。

- `force_field.family`：GAFFまたはGAFF2。
- `resp`：QM backend、method、basis、QM geometry file、ESP条件、RESP fit条件。`source_file` は**ASE1分子の原子順序で9個**の電荷（空白区切り、総和0）。stage入力やQM/ESP原本も保存して出典を追跡します。
- `structure.parameterization_file`：buildが出力した `ethanol_single.sdf` のパス。config相対で `../../runs/ethanol27_initial/ethanol_single.sdf`。
- `structure.ase_atom_order`：GROMACS単分子topologyの原子iが、ASEの何番（0始まり）かを並べた9整数。元素だけで同じ水素同士の対応は確定できません。結合・原子名・QM原子順序を確認します。
- `nonbonded`：電荷・LJモデルに整合する長距離処理、cutoff、分散補正、採用理由。本adapterは明示的なPME + Cut-offを扱います。小さい1.38 nm boxではrlist/rcoulomb/rvdwは0.69 nm未満。この小さいcutoffを物性評価用の推奨値とは扱いません。
- `review`：力場パラメータ/1–4設定の確認と許容最大force。閾値は科学判断として明示し、初期構造を自動で最小化しません。

外部でRESPをフィットした後、既存WorkbenchのAmberTools/ACPYPE pipelineを利用できます：

```bash
"$MM" run -n gmx-md python -m mdworkbench.ethanol_benchmark parameterize --config systems/ethanol_27/system.yaml --output runs/ethanol27_parameters
```

`resp.atom_order_confirmed: true` はQM→電荷→SDFの順序を確認した後だけ設定します。原子順序を変更した場合は電荷も対応させます。`final.mol2`、`parameters.frcmod`、AMBER topology、ACPYPE出力とmetadataを点検してください。RESPフィッティングを未指定条件で自動実施しません。既存の単分子topologyを持ち込む場合も、電荷と原子mappingを検証します。

## 入力が揃った後の操作

`--single-top` にはACPYPEの**1分子**GROMACS `.top` を指定します。下のパスは出力を確認して置換してください。includeが別directoryなら `--include-dir` を追加できます。

```bash
"$MM" run -n gmx-md python -m mdworkbench.ethanol_benchmark prepare --single-top runs/ethanol27_parameters/system.acpype/system_GMX.top
"$MM" run -n gmx-md python -m mdworkbench.ethanol_benchmark preflight
"$MM" run -n gmx-md python -m mdworkbench.ethanol_benchmark run
"$MM" run -n gmx-md python -m mdworkbench.ethanol_benchmark analyze
```

prepareは27分子へのtopology展開、全includeのsnapshot、RESP照合、結合mapping・質量・パラメータの検査、高精度GRO座標/速度、専用MDP、hash、parameter reportを生成します。GROの桁幅を広げ、nm・nm/psへ変換します。`grompp` に `-maxwarn` を渡しません。TPRとprocessed topologyはGROMACSが解釈した入力です。

preflightは元の初期座標で**0ステップのエネルギー/力評価**を行います。非有限値、巨大force、分子間の極端な接近を検出したら停止します。forceの閾値を超えることは、元配置では安全に開始できない可能性を意味します。閾値を引き上げる前に原因を調べます。EMを勝手に実施しません。preflight後も元の速度・座標から1本のNVTを開始し、preflight checkpointから継続しません。

runは0〜600の連続計算です。`nvt.xtc` は5stepごと、quantizationは0.00001 nm。`nvt.edr`→`thermodynamics.xvg` は温度とPotential（XVGのlegendを読み、列順を推測しません）。解析はstep205〜600の80frameだけを許可します。

## 温度の自由度とthermostatの対応

初期system COM運動量を1回除去して、物理的な内部自由度 `3×243−3=726` で瞬間温度を300 Kへ1回だけスケールします。分子ごとのCOM/回転は除去しません。ASEの拘束なし `get_temperature()` は729自由度なので298.765 Kと表示します。

GROMACSはMD中にCOMを繰り返し除去しないよう `comm-mode=none` とします。この場合のcoupling自由度も729であり、GROMACS初期表示温度は約298.765 Kです。初期条件とエンジン内の温度定義を完全に同時一致させていません。GROMACS topology質量とASE質量の小さい丸め差も残り得ます。reportとLOGの自由度・初期温度を確認します。元MACE入力がないため、その速度との一致は未検証です。

ASEのNHC質量が `Q1=ndof kBT tdamp²`、GROMACSの時間パラメータが振動周期を表すことから、同じ自由度を想定すると `tau_t=2π tdamp=0.3141592654 ps` を候補として採用します。これは積分分割や内部状態まで同じという保証ではありません。参考：[ASE NHC source](https://docs.ase-lib.org/_modules/ase/md/nose_hoover_chain.html)、[GROMACS thermostat理論](https://manual.gromacs.org/2022/reference-manual/algorithms/molecular-dynamics.html)、[MDP options](https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html)。

## RDFの共通定義

単位はÅ。分子ごとにOとOHのHを1個ずつ選び、O_i→O_j、O_i→H_jでi≠jだけを数えます。OOは方向付きpairで、各中心の26partnerを含む702pairが1frameの候補です。

`g(r)=全frameのbin pair数 / Σframe[27 × (26/V) × (4π/3)(r_outer³−r_inner³)]`

45bins、r_max=6.762 Å、minimum image。27/Vで規格化した値はこの定義の26/27倍になります。理想的な一様配置を十分平均すると基準は約1です。今回の短い格子由来の系列が1に近づく保証はありません。OH RDFは別分子のOH hydrogenとの接近分布で、角度条件を含む水素結合数ではありません。

MACE trajectoryが将来用意できた場合：

```bash
"$MM" run -n gmx-md python -m mdworkbench.ethanol_benchmark analyze --mace /path/to/reference.extxyz
```

ASEの原子順・分子順を保ち、各frameの`step`情報が必要です。ない場合はstep配列`.npy`を `--mace-steps` で指定します。保存間隔を推測しません。GROMACS/MACE readerは同じRDF関数へÅ座標を渡します。比較元の配列・metadataも検証してください。異なるpotentialの絶対Potential Energyを直接一致させる指標にはしません。

## 考察課題

1. 同じ初期座標・速度でも加速度が異なると、trajectoryはどう変わるか。
2. tdampとtau_tを対応させてもNVT dynamicsが完全一致しない理由は何か。
3. constraints=noneで0.5 fsを採用する理由を結合振動から説明する。
4. 200 fsのproductionと80枚の隣接frameを、80個の独立標本と考えてよいか。
5. OOとOH RDFのピークは何の距離を示すか。初期格子のピークはどこにあるか。
6. 26/Vと27/Vで規格化がどれだけ変わるか、コードで確認する。
7. 条件を揃えた比較で、依然として何が異なるか（mass、DoF、precision、thermostat、長距離処理）。
8. 平衡構造・実験RDF・拡散係数を主張するため、何を追加する必要があるか。

全300 fs（equilibration100 fs、production200 fs）は操作とモデル差の教材です。平衡液体構造、収束RDF、信頼できる物性を得たとは主張しません。
