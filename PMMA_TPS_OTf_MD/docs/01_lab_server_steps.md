# 10月29日：研究室サーバーでの開始手順

この文書のコマンドは**研究室Linuxサーバー用**です。Macでg16/実分子RESP/MDを実行する指示ではありません。以降の`python`は、サーバーで確認したPython環境のコマンドに読み替えます。

## 1．フォルダを置き、実行環境を調べる

`PMMA_TPS_OTf_MD/`を丸ごと、研究室のwork等へコピーします。Pythonコードは自分のフォルダからの相対位置でファイルを解決するため、Macの絶対パスを変更する必要はありません。元のv2仕様書は参考用として一緒に持っていくと便利ですが、計算スクリプトの実行にrepo全体は不要です。

```bash
cd PMMA_TPS_OTf_MD
python scripts/server.py probe
```

`runs/server_probe.json`で、g16/formchk、AmberToolsのespgen/resp/antechamber/parmchk2/tleap/sander、ACPYPE、GROMACS、qsub、Pythonパッケージを確認します。未導入なら勝手に別法へ変更せず、不足を記録して環境を整えます。既に再現に成功した環境を優先してください。

`probe`が示すCPU数はマシン全体の値であり、その数を自分のジョブで使える意味ではありません。

## 2．Gaussian用の資源設定を確定する

```bash
cp config/server.example.json server.local.json
```

`server.local.json`を編集し、次を確認した値で埋めます。

| キー | 設定するもの |
|---|---|
| `configured` | 全項目を確認した後に`true` |
| `nproc` | 1ジョブに許可されたCPUコア数。Gaussianの`%nprocshared`にも使う |
| `memory_gb` | Gaussianの`%mem`。ジョブ全体のメモリ上限ではない |
| `pbs_memory_gb` | PBSが確保する全メモリ。`memory_gb`より大きくし、OS等の余裕を設ける |
| `queue` / `walltime` | 利用可能なPBS queueと`HH:MM:SS` |
| `gaussian_setup` | g16環境を有効にするshellファイル。既に適切なら`null`可 |
| `amber_setup` | AmberTools/Python等を有効にするshellファイル。同様に`null`可 |
| `python_command` | NumPy/RDKit/ParmEd/pyedr等が使えるPythonの実体 |
| `g16_command` | g16の起動コマンドまたは実行ファイル |
| `scratch_root` | 計算ノードから書き込める、容量のある既存の絶対パス |
| `esp_method` | 今回は`"mk"` |

環境を`module load`等で有効にする研究室なら、それを記載した自分用のshellファイルを作って`gaussian_setup`等へ指定します。元の成功済みPBS例と、queue・`select`・CPU/memory記法を照合してください。テンプレートはPBS Pro系の`select=...`記法です。研究室側の記法が異なれば、実際のジョブ例に合わせて修正が必要です。

## 3．実行用Gaussian入力とPBSジョブを展開する

```bash
python scripts/server.py render --profile server.local.json --esp-method mk
```

`runs/qm/<job_id>/`に次が生成されます。まだジョブ投入はしません。

```text
opt.gjf       # 初期座標からRHF/6-31G*構造最適化
esp.gjf       # optimized.chkを読み、同じ構造でMK ESP
job.pbs       # opt→正常終了確認→ESP→formchk→espgen→整合検査
job.json      # 分子・配座・用途・設定の記録
```

同じ`runs/qm/`があれば再生成を停止します。入力変更のために実行途中の出力を上書きせず、別のフォルダコピー/runを作るか、未投入であることを確認した上で整理してください。

## 4．まずOTf⁻だけを投入する

```bash
qsub runs/qm/OTF_c01/job.pbs
```

投入時に表示されたjob IDを記録し、`qstat`で状態を確認します。`opt.log`が存在するだけで完了とは判断しません。

ジョブが正常に終わると、次が揃います。

- `opt.log`：最適化の収束、SCF収束、Normal termination。
- `esp.log`：MK ESP計算のNormal termination、ESP fit中心・点。
- `optimized.chk` / `optimized.fchk`：同じ構造・波動関数。
- `qm.esp`：Amber RESP形式の核座標・格子座標・静電ポテンシャル。
- `optimized.sdf`：原子順を維持した最適化後構造。
- `result.json`：原子数・電荷・多重度・構造対応・ESP単位/向き・hashの検査結果。

`result.json`がない場合は、最初に失敗したログを読んで原因を解決します。エラー後に次の重いジョブをまとめて投入しないでください。失敗したジョブのscratchは調査用に残り、成功時だけそのジョブ用scratchを削除します。

OTfの実測時間・メモリ・研究室の同時実行制限を確認して、TPS_c01〜c03、PMMA3_c01〜c03、PMMA5_c01〜c03へ進みます。PMMA19の全鎖QMジョブはありません。

## 5．RESPを行う

以下の計算コマンドは、研究室の利用規則に従って計算ノードのPBSジョブまたは対話割当内で実行します。PBSの投入方法は元の成功例を利用してください。

OTfを先行して確認できます。

```bash
python scripts/fit_resp.py fit --group OTF
```

TPSの3配座、PMMA3/5の各3配座が揃ったら：

```bash
python scripts/fit_resp.py fit --group TPS
python scripts/fit_resp.py fit --group PMMA
```

処理は実ESPを検査し、点数に応じた重みを生成し、AmberToolsの`resp`を呼びます。PMMAではtemplate転写モデルと、比較用の分子ごとの非転写RESPを評価します。fitに使わないc03の誤差も計算します。

結果：

```text
runs/resp/OTF/     # stage1.in/out/qout、ESP、実行記録、review.json
runs/resp/TPS/     # 同上
runs/resp/PMMA/    # stage1/2、比較用fit、転写・評価のreview.json
runs/resp/charges/ # OTF.charges、TPS.charges、PMMA19.charges
```

`review.json`の`validation`、電荷の和、部分電荷の分布、評価配座のESP誤差、転写による誤差増加を確認します。3つのPMMA配座が最適化後に重複した場合は停止します。異なる初期配座を作って再計算し、評価用データを確保します。

実データのESP再現を図で確認するには、計算済みの群について次を実行します。

```bash
python scripts/plot_resp.py OTF
python scripts/plot_resp.py TPS
python scripts/plot_resp.py PMMA
```

各`runs/resp/<群>/figures/`へ、QM ESP対RESP ESPの散布図と残差分布を保存します。評価用c03も確認し、fit用だけが良く合う結果を見逃さないようにします。

確認が済んだ群だけ、その`review.json`の`approved`を`true`にします。これは実データを見て採用するためのゲートで、未計算のファイルへ設定するものではありません。

### RESPの出力桁で止まった場合

採用するRESP版によっては`qout`が小数6桁です。10桁で保存し直しても元の精度は増えません。19量体へ繰り返し転写した結果、丸め誤差が許容値を超えたら処理を停止します。一律オフセットで中性にせず、RESPの高精度出力方法、拘束付き再fit、丸めの影響を検討してください。この問題が生じた段階の入力・出力を保存し、次のMDへ進めないでください。

## 6．GAFF2のトポロジーを生成する

```bash
python scripts/forcefield.py build OTF
python scripts/forcefield.py build TPS
python scripts/forcefield.py build PMMA19
```

実RESP電荷を読み込んだmol2、frcmod、Amberのprmtop/inpcrd、ACPYPE変換後のGROMACS topologyが生成されます。`runs/forcefield/<分子>/`へ保存します。未生成の電荷・未承認fit・未審査の欠損パラメータでは止まります。

`build_report.json`、frcmod、原子型、S周辺、主鎖二面角、improper、1–4を検査します。問題を解決したら、各`build_report.json`の`parameter_review_approved`を`true`にします。ゼロ項や類推値の意味を調べずに採用しないでください。

## 7．AMBER/GROMACSの成分別エネルギーを比べる

まず単分子を確認します。`--gmx`には実機で確認したコマンド名または実体を指定できます。

```bash
python scripts/check_energy.py OTF --gmx gmx
python scripts/check_energy.py TPS --gmx gmx
python scripts/check_energy.py PMMA19 --gmx gmx
```

両engineに**同じ精度の座標**を渡し、bond/angle/dihedral/LJ-14/Coulomb-14/LJ/Coulombを比較します。これは0-stepの変換検査です。MDの物性計算ではありません。

次に、主鎖二面角を変えたPMMA19の2構造と、TPS–OTfの距離・向きを変えた6構造を準備します。

```bash
python scripts/conversion_cases.py prepare
```

`runs/conversion_cases/matrix.json`の`cases`に列挙された8フォルダについて、計算ノード上で次を実行します。`CASE_DIRECTORY`は実際の1フォルダへ置き換えます。

```bash
python scripts/check_energy.py --case CASE_DIRECTORY --gmx gmx
```

8件が完了したら：

```bash
python scripts/conversion_cases.py summarize
```

全成分・全ケースとパラメータレビューに根拠が揃ったときだけ、`reports/conversion_review.json`が合格になります。差が出た場合は、力場項・単位・1–4・LJ混合・座標精度などを調べます。合計値の相殺で合格にしません。

## 8．混合系の入力へつなぐ

```bash
python scripts/forcefield.py assemble
python scripts/pack_box.py --gmx gmx
```

`runs/adopted_input/`へ実際の`.itp`・`.top`・単分子`.gro`を統合し、`runs/initial_box/init.gro`を作ります。PMMA100→TPS83→OTf83の順に挿入し、各段の実挿入数と最終32,186原子を確認します。

次にEMの入力とコマンドを準備します。

```bash
python scripts/prepare_md_stage.py opt --gmx gmx
```

`runs/md/opt/commands.sh`は**コマンドの準備まで**で、自動実行・投入しません。元の成功済みPBSジョブを基に、同じGROMACS/MPI環境と資源を確保し、計算ノードで実行します。事前の0-step力評価、Fmax・最小距離・grompp警告等の確認も必要です。

0-step用の`mdp/preflight.mdp`も準備しています。出力先を作り、計算ノード上で例えば次を実行します。

```bash
mkdir -p runs/preflight
gmx grompp -f mdp/preflight.mdp -c runs/initial_box/init.gro -p runs/initial_box/topol.top -o runs/preflight/preflight.tpr -po runs/preflight/effective.mdp -pp runs/preflight/processed.top
gmx mdrun -deffnm runs/preflight/preflight
```

これは時間を進める本計算ではなく、初期構造の力・エネルギーの検査です。生成したTRRの力の単位、最大力の原子、最短距離を確認し、NaN/Infや異常な近接があれば構造・力場へ戻ります。問題がないことを確認してからoptへ進みます。

optの検証後は、実際の出力を指定してheatを準備します。

```bash
python scripts/prepare_md_stage.py heat --parent-gro runs/md/opt/opt.gro --gmx gmx
```

その後は各段の実GRO/CPTを渡します。`.partXXXX`付き出力があるときは、実際のファイル名を指定してください。

```bash
python scripts/prepare_md_stage.py eq_1 --parent-gro ACTUAL_HEAT_GRO --parent-cpt ACTUAL_HEAT_CPT --gmx gmx
```

予定とcheckpointのstepが違えば停止します。上流の延長後に古い`init-step`や冷却時刻をそのまま使わず、スケジュールを再生成・検証します。

## 明日の最初の到達点

最初から全77.2 nsを完了させるのではなく、まず次を確認してください。

1. 実行環境と資源が確認できた。
2. OTfの最適化とMK ESPが正常終了した。
3. 核座標・原子順・電荷・多重度・ESP点数の検査に合格した。
4. OTfのRESPを実行し、部分電荷と和−1 eを説明できた。
5. 残りのQMジョブの資源・計算時間を実測に基づいて決めた。

ここまでのログが揃えば、その後のPMMA転写・力場検証へ、実結果を見ながら進められます。
