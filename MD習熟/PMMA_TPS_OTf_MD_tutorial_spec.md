# PMMA / TPS⁺ / CF₃SO₃⁻ レジストモデル MD 習熟チュートリアル仕様書

版：1.0 / 作成日：2026-10-27（JST） / 状態：仕様策定済み・新規系の構築と計算は未実施

## 1. 目的と、この仕様書の使い方

参照 [MD_output_guide.ipynb](MD_output_guide.ipynb) の結果を再現できた利用者が、**分子を自分で定義し、力場・初期配置・MD・冷却・解析を一から構築する**ための実装仕様である。実行先は参照計算を再現した環境を使う。

到達点は、TPSカチオン・PMMA・トリフラートアニオンからなる中性のバルクモデルについて、次を自分で説明し、同じ入力から再実行できることである。

1. 構造と力場を準備し、`.gro` / `.itp` / `.top` / `.mdp` の対応を確認する。
2. `gmx insert-molecules` で初期配置を作り、最小化とNPTによる密度調整を行う。
3. 高温MDの後、冷却・室温NPT再平衡化・室温Production MDを行う。
4. 密度、エネルギー、MSD、拡散係数の評価可否、RDFを検証する。
5. イオン対、クラスター、イオンとポリマーの接触、ポリマー鎖・局所運動を解析する。
6. 数値計算の完走、観測量の安定、統計の十分さ、実材料への妥当性を区別する。

**「必ず実験物性が再現される」という保証はしない。** 習熟を確実に進めるため、各Phaseに入力・成果物・合格基準・失敗時の分岐を設ける。未収束や拡散未観測も、根拠とともに判定できれば学習成果になる。数値結果を得られない段階を、得られたことにして先へ進まない。

本書の数値は、参照Notebookの保存情報、そこからの計算、本教材で採用する設計値に区別する。**本書内のスクリプト名・Notebook名・設定形式は、特記以外は今後実装するもの**であり、現在使えるコマンドと誤認しない。

## 2. 参照から確認できたことと、実行環境の扱い

### 2.1 参照情報

次は参照Notebookのセルと保存出力を読んだ結果であり、元のバイナリ・入力ファイルをこのMacで再検査した結果ではない。利用者の「再現できた」という申告を前提に、再現作業のやり直しは要求しない。

| 項目 | 参照Notebookで確認した値 | 新規系への使い方 |
|---|---|---|
| 元ポリマー | `POL.000` / 304原子 / 約1,952.39 g/mol | PMMAの**1鎖の分子量**を合わせる目安 |
| 元ポリマーの組成表示 | `C106H166O32` | PMMAそのものではない。分子構造を流用しない |
| 元カチオン | TPS / `C18H15S` / 34原子 / +1 e | 同じ化学種を新規作成する |
| 元アニオン | F1 / `C14H20FO6S` / 42原子 / −1 e | CF₃SO₃⁻とは別物。F1の構造・電荷を流用しない |
| 分子数 | ポリマー100鎖、TPS83個、F1 83個 | 基本案は鎖数とイオン対数を維持 |
| 初期箱 | 10 × 10 × 10 nm³ | 基本案として維持。挿入可否は新規系で確認 |
| 最小化 | `steep`、最大50,000 step、`emtol=100` | 終了理由とFmaxを読む練習に使う |
| NPT | 1400 K、100→10→1 bar、0.5→0.5→5 ns | 段階構成を継承する |
| 高温本計算 | 1400 K、1 bar、20 ns、1 fs/step | 高温での混合・解析練習として継承する |
| 制御 | V-rescale / C-rescale、τT=1 ps、τP=5 ps | 基本設定として継承する |
| 元の段間接続 | 各段で`gen_vel=yes`、`continuation=no` | 新規系は最初だけ速度を生成し、その後は継続する |
| 先輩の計算 | GROMACS 2022.5、16 MPI rankというログ | 過去の実行条件として保存 |
| 利用者の再現計算 | GROMACS 2026.3という記載・保存出力 | **新規系の実行環境の第一候補** |
| 実行方式 | PBS、16 CPU要求、利用者側は30時間walltime追記との記載 | 現在の再現成功ジョブを実機で読み直して確定 |

参照の高温密度約0.447 g/cm³やRDFピークは、組成と力場が異なる新規系の合格値には使わない。Notebook内の「平均密度差0.01 g/cm³なら同じ」という目安も、そのまま一般化しない。

### 2.2 「同じ環境」の定義

この作業ディレクトリはMac上にあるが、参照Notebookの保存出力にはLinux/HPC上の実行パスとPBSジョブが記録されている。Mac側の既存GROMACS 2025.4環境を、再現成功環境と同一視しない。

Phase 0で、利用者が再現に使った環境から次を採取し、`environment_report.json` と `environment_report.txt` に保存する。

- OS、CPUアーキテクチャ、メモリ、空き容量、作業用領域、保存先。
- GROMACS実行ファイルの実体と `gmx --version`、MPI/thread-MPI、精度、SIMD、GPU対応。
- Python interpreter、Jupyter kernelspecの`argv[0]`、パッケージ一覧。
- AmberToolsの`antechamber` / `sqm` / `parmchk2` / `tleap`、ACPYPEの実体とバージョン。
- RDKit、NumPy、pandas、Matplotlib、MDAnalysis、pyedr、nbformat、PyYAMLの利用可否。
- 成功したPBSスクリプト、queue、resource要求、launcher、環境変数、PATH、実際のrank/thread数。
- 実行バイナリとモジュール一覧。バージョン文字列だけで一致としない。

Notebookのカーネル表示名は `Python (GROMACS MD Workbench)` / IDは `gmx-md`。ただし、同名でもマシンごとに実体が異なるため、**実行先で `sys.executable` とtool実体を照合**する。既存の環境検査コードは再利用できるが、既存`probe_server.sh`はSlurm中心なのでPBSの`qstat`や成功済みジョブの確認を追加する。

まず既存環境で実行可能かを検査する。不足ツールがあれば、追加内容・保存先・容量・既存環境への影響を提示してから導入する。今回の仕様書作成だけでは、インストール、環境更新、HPCへの接続・ジョブ投入を行わない。

## 3. モデルの定義

### 3.1 基本組成

本書ではTPSを**トリフェニルスルホニウム（triphenylsulfonium）カチオン**と解釈する。略称だけでなく、構造式・形式電荷・分子式で固定する。

| 化学種 | トポロジー分子名 | 残基名の案 | 電荷 | 基本個数 | 構造の定義 |
|---|---|---|---:|---:|---|
| PMMA19 | `PMMA19` | `PMM` | 0 | 100鎖 | MMAをhead-to-tail結合した19量体、両末端H終端 |
| TPS⁺ | `TPS` | `TPS` | +1 | 83個 | Sにフェニル基3つが結合したカチオン |
| CF₃SO₃⁻ | `OTF` | `OTF` | −1 | 83個 | トリフラート。トリフルオロメタンスルホン酸イオン |

構造生成の入力例：TPSは `c1ccc([S+](c2ccccc2)c2ccccc2)cc1`、OTfは `O=S(=O)([O-])C(F)(F)F`。SMILES文字列は分子グラフ検査を通してから採用する。CF₃SO₃Hや、中性のS原子を持つ別のTPSを作らない。

初期箱は立方体10 nm角、PBCはxyz、溶媒・水・基板・真空層は入れない。1鎖の分子量と箱全体の総分子量は区別する。

### 3.2 PMMAの重合度・末端・立体規則性

繰り返し単位は `–CH₂–C(CH₃)(C(=O)OCH₃)–`、組成はC₅H₈O₂、約100.117 g/mol。モノマーのC=Cを残したまま単に19個並べた構造にしない。

採用案は**重合度19、直鎖、head-to-tail、両端の未結合価をHで終端**する。したがって組成はC₉₅H₁₅₄O₃₈、全水素を含め287原子、分子量は約1,904.24 g/molとなる。参照1,952.39 g/molとの差は約−2.47%。20量体は約2,004.36 g/mol、差約+2.66%なので、H終端の条件では19量体を選ぶ。

これは分子量を合わせた**単分散オリゴマー教材**であり、長鎖の実用PMMAの分子量分布・絡み合い・Tgをそのまま表さない。

立体配置は、seed 42で生成するランダムな固定立体配列を基本とする。実際に立体中心となる原子を分子グラフから判定し、atom-map番号、絶対配置、meso/racemo dyadを保存する。「R/Sをランダムに指定した」という記録だけでアタクチック性を保証しない。鎖末端などで立体中心でない原子も区別する。

最初の教材では、この**1つの固定立体配列のトポロジーを100鎖に使用**する。これはアタクチック鎖の一実現を複製した近似である。立体配列を複数に増やす検討は発展課題とし、基本系の途中で勝手に変更しない。各鎖の向き・位置は挿入時にランダム化する。初期配座が共通である影響は、高温でのRg・二面角・接触の時間変化で確認する。

### 3.3 分子数比と質量分率

参照と同じ個数比を採用すると、次の設計値になる。値は一般的な原子量からの概算で、最終値は**採用したトポロジーの質量**から再計算する。

| 量 | 基本案の値 |
|---|---:|
| TPS⁺の分子量 | 約263.38 g/mol |
| OTf⁻の分子量 | 約149.06 g/mol |
| 1塩対の式量 | 約412.44 g/mol |
| 総原子数 | **32,186** |
| 分子数 | 266（100鎖＋83＋83） |
| 総電荷 | 0 e |
| イオン対数 / PMMA鎖数 | 0.83 |
| イオン対数 / MMA単位数 | 約0.04368 |
| 塩の質量 / 全質量 | **約15.24 wt%** |
| 塩の質量 / ポリマー質量 | 約17.98 wt% |
| 初期密度（10 nm角） | 約0.3731 g/cm³ |

塩の比率を実材料に合わせる依頼がないため、今回は個数比を維持する。「参照と同じ個数比」と「参照と同じ塩の重量比」は両立しないことを説明する。別の配合を試すときは別config・別runとして扱う。

### 3.4 物理モデルの適用範囲

本教材の標準は、固定電荷・固定結合の非反応性全原子MDとする。

- 観測するのは、未露光バルク中のPAGイオンの配置、会合、移動、ポリマーとの局所接触と鎖運動である。
- TPSの光分解、酸の発生、プロトン移動、酸触媒脱保護、架橋、溶解、パターン形成は計算しない。
- **OTf⁻の移動を「光酸の拡散係数」と呼ばない。** 酸を生成する反応もプロトンもモデルにない。
- PMMA単独は、参照の酸脱保護型ArF樹脂と同一材料ではない。これは指定3成分によるレジスト関連モデルである。
- 1400 Kは参照に合わせた計算上の高温処理である。実材料の熱履歴ではなく、力場の通常の較正範囲からも外れ得る。結合が切れないことは、高温の物性や配置が正しいことを意味しない。
- 分解・分極・電荷移動を無視する誤差、末端・有限サイズ・冷却速度・力場による誤差を最終報告に残す。

1400 Kで数値不安定や極端な分子内変形が起きたら、時間刻みと力場を検査する。無断で温度を下げて「同条件」とはしない。600–800 K程度から開始する代替プロトコルは、高温のモデル外挿を減らせる一方、緩和が遅くなり参照から条件が変わるため、別の設計判断として提示する。

## 4. 実装成果物と保存構成

既存エタノール教材と参照Notebookは保存し、新規系専用の構成を追加する。共有ライブラリを変更する場合も既存RESP入力の必須条件を緩めない。

```text
examples/pmma_tps_otf/
  README.md                         # 実行先・カーネル・最短の開始手順
  PROGRESS.md                       # Phase状態、証拠、次の1手
  decisions.md                      # 採用条件・変更理由・未解決事項
  config/system.yaml                # 物理モデルと組成の正本
  config/protocol.yaml              # ステージ・時間・保存間隔の正本
  notebooks/
    00_environment_and_plan.ipynb
    01_molecular_structures.ipynb
    02_force_field.ipynb
    03_initial_packing.ipynb
    04_minimization_and_preflight.ipynb
    05_high_temperature_md.ipynb
    06_cooling_and_room_temperature.ipynb
    07_standard_analysis.ipynb
    08_resist_molecular_behaviors.ipynb
    09_validation_and_report.ipynb
  scripts/                          # Notebookからも呼ぶ処理の実装
  templates/pbs/                    # 実機確認後のPBSスクリプト
  templates/report.md
systems/pmma_tps_otf/
  structures/                       # canonical SDF、atom_map、site_map、立体配置
  parameters/                       # mol2、frcmod、prmtop、採用gro/itpと検査記録
  topology/                         # 共通defaults/atomtypes、topol.top
  mdp/                              # configから生成した全段の入力
runs/pmma_tps_otf/<run_id>/
  manifest.json                     # 入出力実パス・hash・実行条件・親run
  environment/                      # 実機検査・package一覧
  00_pack/ 01_opt/ 02_heat/ 03_eq_1/ ...
  analysis/                         # CSV、XVG、図、解析条件
  executed_notebooks/               # 実行済みNotebook。元教材と分離
  reports/                          # 合格判定、性能、最終報告
```

`runs/`は既存のignore方針を踏襲する。HPCでの大容量出力は、その環境のwork/scratchへ保存し、configとmanifestで対応付ける。scratchの保持期限を調べ、削除前に最終構造・checkpoint・入力・解析結果・必要な軌跡を保存先へ退避する。

各stageは独立ディレクトリで実行する。`mdout.mdp`の上書きを避け、`<stage>_effective.mdp`、`<stage>_processed.top`、gromppのstdout/stderrを残す。原本の構造・力場は採用後にhashで固定する。

### 4.1 設定の正本と最小例

実装時に`system.yaml`へ最低限次を保持する。これは**新規に実装する設定形式の案**であり、既存CLIへそのまま渡せる設定ではない。計算の前にschema検査を行い、未設定の実行先・tool・力場レビューを暗黙の既定値で補わない。

```yaml
schema_version: 1
system_id: pmma19_tps_otf
model:
  polymer:
    degree: 19
    linkage: head_to_tail
    end_caps: [H, H]
    stereochemistry: seeded_fixed_sequence
    stereo_seed: 42
    conformer_seed: 43
    count: 100
  cation: {identity: triphenylsulfonium, charge_e: 1, count: 83}
  anion: {identity: triflate, charge_e: -1, count: 83}
  initial_box_nm: [10.0, 10.0, 10.0]
  pbc: xyz
force_field:
  family: gaff2
  charge_model: am1_bcc
  charge_scale: 1.0
  review_status: pending
execution:
  target: reproduction_success_environment
  scheduler: pbs
  run_root: null
  gmx_bin: null
  python_bin: null
  queue: null
  launcher: null
```

`protocol.yaml`にはstageごとのT/P/ensemble、`duration_ps`、dt、出力間隔、冷却の相対点列、延長履歴を持たせる。MDP・step数・累積時刻はそこから生成する。固定体積stageには`pressure_bar: null`とし、300 KのNVTへ「1 barを制御した」という情報を付けない。未知の`null`は実行前に解決し、力場レビュー未完了なら本番投入を禁止する。

## 5. Phase 0：環境確認と計画の固定

**入力：** 本仕様、参照Notebook、再現成功環境とそのジョブファイル。

**実装内容：** 環境検査、参照条件の抽出、基本config、ディレクトリ、進捗台帳を作る。Notebookの保存出力にある絶対パスを新規コードへ直書きせず、`repo_root` / `run_root` / `reference_root` / `gmx_bin`で管理する。

**成果物：** `00_environment_and_plan.ipynb`、環境レポート、config、`PROGRESS.md`、`decisions.md`。

**合格基準：** 実行先のPython・GROMACS・PBSの対応が確認でき、依存不足と入出力先が明記されている。構造・力場・長時間計算の未実施状態が区別されている。

HPCにCodexから接続できない場合は、実行コマンドとPBSスクリプトを用意し、利用者が実行したログを検査して進める。ローカルの成功をHPCでの成功として記録しない。

## 6. Phase 1：分子構造を一から作る

### 6.1 生成手順

1. RDKitでTPSとOTfの結合・形式電荷を指定し、明示的H付きSDFを保存する。
2. PMMAの繰り返し単位を、結合点を明示した分子グラフ操作で19個結合する。末端H、主鎖の飽和、エステルのC=O、総形式電荷0を検査する。
3. PMMAの立体配置を固定し、seed 43で初期3D配座を生成する。配座生成・分子内最適化のmethod、収束、失敗も保存する。
4. canonical SDFを正本とし、mol2/gro等へ変換する。全形式間のatom-mapを作る。
5. 3D表示と、分子式・質量・電荷・結合数・結合次数・立体配置の表を作る。

RDKitのMMFF等による構造調整は、初期配座作成の補助である。TPSなどでパラメータ非対応なら、既定法の成功を装わず、後述のGAFF2最小化等へ明示的に切り替える。最終MD力場と同じものを使ったと誤記しない。

### 6.2 原子と解析サイトの対応

`atom_map.csv`に化学種、canonical atom ID、元素、形式電荷、結合相手、mol2 ID、GROMACS原子番号・原子名を保存する。元素を原子名の先頭1文字だけから推定しない。

`site_map.json`には少なくとも以下を定義する。

- TPSの中心S、芳香環C/H。
- OTfのS、3つのスルホン酸O、CF₃のC/F。
- PMMAの主鎖C、カルボニルO、エステル単結合O、側鎖メチル、末端、繰り返し単位ID。
- 分子IDと鎖ID。**残基と分子は別概念**として扱う。

`.gro`の番号は桁制限で巡回する場合がある。分子数を`residの一意数`だけで判定せず、トポロジーの分子境界とatom-mapで検査する。

**成果物：** `01_molecular_structures.ipynb`、3種のcanonical SDF、立体配置記録、map、構造図。

**合格基準：** PMMA19=287原子、TPS=34原子、OTf=8原子、形式電荷0/+1/−1、意図した結合グラフ、立体配置、変換前後の原子対応が一致する。PMMA中に未反応モノマーのC=Cがない。

**失敗時：** 分子量を合わせるために原子を足し引きしない。結合点・末端・Hの付加を修正して再生成する。

## 7. Phase 2：力場・電荷・GROMACS入力を作る

### 7.1 標準経路と変更条件

**本教材の基本案はGAFF2＋AM1-BCC、非分極・全電荷（0/+1/−1）モデル**とする。AmberToolsで3種を同じ規則により処理し、LEaP→ACPYPEでGROMACSへ変換する。これは実行負担を抑えて系の構築全体を学ぶための設計選択で、イオン性PMMA系に対する精度保証ではない。

AM1-BCCとRESPは別の電荷モデルである。RESPを使う場合は、構造最適化条件、QM法・基底、ESP格子、配座、電荷制約、fit条件を記録して、**全3種を含む別の力場variant**として扱う。RESPの名称をAM1-BCC出力へ付けない。GAFFの背景と電荷モデルの関係は[Amber公式GAFF解説](https://ambermd.org/antechamber/gaff.html)、実行オプションは[Antechamber公式解説](https://ambermd.org/antechamber/ac.html)と実機の`-h`で確認する。

既存`mdworkbench.parameterization.amber`はRESP入力を要求する設計である。本教材のAM1-BCC経路は専用モジュール・明示的な設定として追加し、既存エタノール系で不足RESPを自動BCC補完する変更はしない。

PMMA19は、末端を含む**1本の鎖全体に対してAM1-BCCを計算する**。単量体の電荷を19回複製して鎖を作らない。鎖サイズでSQM等が不収束・資源超過になった場合は、原因と費用を報告し、検証済みフラグメント法または文献パラメータへの変更を提案する。無断で電荷ゼロにしない。

### 7.2 変換パイプライン

実装はcanonical SDFから、次の生成と記録を行う。

```text
SDF（形式電荷・原子対応を固定）
  → antechamber（GAFF2型、AM1-BCC、nc=0/+1/−1、multiplicity=1）
  → charged mol2 + SQM/antechamberログ
  → parmchk2（GAFF2）→ frcmod
  → tleap（leaprc.gaff2、check、prmtop/inpcrd）
  → ACPYPE（既存Amber topologyを変換し、電荷を再計算しない）
  → 単分子gro/itp → 共通力場 + topol.top
```

各コマンドの引数、終了コード、stdout/stderr、使用データファイルのhash、tool versionを保存する。自動推定されたパラメータは、生成されたことと妥当性を検証したことを分ける。

### 7.3 必須レビュー

- SQMと電荷計算の収束、部分電荷の一覧、原子対応、分子電荷の合計を確認する。
- TPS中心Sの原子型とC–S結合・C–S–C角・二面角、OTfのS–O・S–C・C–Fを重点確認する。元TPSの`s4`等のラベルを正解として転記しない。
- PMMAの主鎖、エステル、末端、立体配置を保つimproperと主鎖二面角を確認する。
- `ATTN`、未定義、ゼロの力定数、要確認の類推パラメータを抽出する。ツールのexit code=0だけで合格にしない。
- 型ごとのσ・ε、質量、全bond/angle/dihedral/improper、exclusion・1–4 pairを検査する。正当に複数項を持つ二面角と誤重複を区別する。
- 不審なSの型や主要な主鎖二面角を解決できなければ停止する。文献値または小分子QMの構造・回転プロファイルで検証する追加計画を作る。
- この段階の単分子最小化と短いMDで、結合・角・立体配置の破綻がないか確認する。これはイオン会合や凝縮相物性の定量精度検証とは区別する。

全原子電荷を8桁程度以上で保持し、分子電荷誤差は目安`1e-5 e以下`、混合系の総電荷誤差は`1e-3 e以下`を要求する。大きな誤差はatom-map・変換・入力電荷を調べる。丸めの補正が必要なら元の高精度電荷から再出力し、無関係な原子に電荷を押し込まない。

### 7.4 トポロジーの統合

```text
forcefield.itp       # [ defaults ] を1回だけ
atomtypes.itp        # 3種で使う全原子型。分子定義より前
PMMA19.itp          # 分子の定義だけ
TPS.itp
OTF.itp
topol.top
```

GAFF2の変換に対応する`nbfunc=1`、`comb-rule=2`、`gen-pairs=yes`、`fudgeLJ=0.5`、`fudgeQQ=0.833333...`を、実際のAmber topology・変換出力と照合する。異なるLJ値を持つ同名atomtypeは上書きせず停止する。同じ値なら重複を整理し、変更記録を残す。

`topol.top`の最小構成は次とする。実装時に生成したファイル名に合わせる。

```text
#include "forcefield.itp"
#include "atomtypes.itp"
#include "PMMA19.itp"
#include "TPS.itp"
#include "OTF.itp"

[ system ]
PMMA19 TPS triflate bulk model

[ molecules ]
PMMA19  100
TPS      83
OTF      83
```

**成果物：** `02_force_field.ipynb`、mol2/frcmod/prmtop/inpcrd、採用gro/itp/top、電荷・項数・型・変換対応のレビュー表、`force_field_review.md`。

**合格基準：** 3分子の電荷・原子対応と全力場項を検査済み、未解決パラメータなし、混合系の`grompp`が未承認警告なしで通る。導入した力場の由来・利用条件を記録する。電荷スケーリングは基本案ではしない。

## 8. Phase 3：初期配置を作る

参照と同じ`gmx insert-molecules`を使う。[公式説明](https://manual.gromacs.org/documentation/2026.3/onlinehelp/gmx-insert-molecules.html)のとおり、これは同じ入力分子を回転・並進して挿入する処理で、鎖の内部配座を探索する処理ではない。

以下は、採用した単分子groを作業ディレクトリへコピーした後の例である。実装では引数リストでコマンドを組み、終了コードを検査する。

```bash
gmx insert-molecules -box 10 10 10 -ci PMMA19.gro -nmol 100 -try 10000 -seed 1 -o init0.gro
gmx insert-molecules -f init0.gro -ci TPS.gro -nmol 83 -try 10000 -seed 2 -o init1.gro
gmx insert-molecules -f init1.gro -ci OTF.gro -nmol 83 -try 10000 -seed 3 -o init.gro
```

各挿入のログ・実挿入数・原子数・分子順序・箱・密度・最短分子間距離を検査する。`-nmol`の要求数が全部入ったことを、出力ファイルの存在だけで判定しない。

挿入半径は力場のLJパラメータと同じ意味ではない。`vdwradii.dat`、atom名からの半径割当、`-radius` / `-scale` / `-rot`の**使用値**を保存する。入らないときは、まず配座・箱サイズ・試行回数を検討し、半径だけを極端に下げて押し込まない。箱を変えた場合は初期密度とconfigを更新する。別packerの採用は変更理由を記録する。

**成果物：** `03_initial_packing.ipynb`、`init0.gro` / `init1.gro` / `init.gro`、3段のログ、検査表、初期配置の3D表示。

**合格基準：** 100/83/83分子、32,186原子、分子順序PMMA→TPS→OTF、総電荷0、意図した箱、明らかな分子内破綻なし。距離検査に加えて次Phaseの0-step力評価を通す。乱数seedの固定だけで異なるGROMACS版間の座標一致は保証しない。

## 9. 全MD段階の計算設計

### 9.1 基本プロトコル

`prod`は参照と対応する**高温NPTのProduction**、`prod_rt`は追加する**300 KのProduction**とする。両方を残す。室温NPTで密度を決めた後、室温の輸送解析は箱変動による扱いを単純にするためNVTを基本にする。

| 段階 | アンサンブル・条件 | 基本時間 | step数（1 fs） | 主な目的 |
|---|---|---:|---:|---|
| `opt` | Steepest Descent、PME | 時間なし | 最大50,000 | 重なりと過大な力を解消 |
| `heat` | NVT、300→1400 K | 0.2 ns | 200,000 | 初回速度生成、段階的加熱 |
| `eq_1` | NPT、1400 K / 100 bar | 0.5 ns | 500,000 | 初期の低密度箱を圧縮 |
| `eq_2` | NPT、1400 K / 10 bar | 0.5 ns | 500,000 | 目標圧力を段階的に下げる |
| `eq_3` | NPT、1400 K / 1 bar | 5 ns | 5,000,000 | 高温で体積・構造の安定を検査 |
| `prod` | NPT、1400 K / 1 bar | 20 ns | 20,000,000 | 高温の構造・運動の参照 |
| `cool` | NPT、1400→300 K / 1 bar | 11 ns | 11,000,000 | 100 K/nsで連続冷却 |
| `eq_rt` | NPT、300 K / 1 bar | 最初10 ns | 10,000,000 | 室温の密度・構造の安定を検査 |
| `nvt_rt` | NVT、300 K | 0.2 ns | 200,000 | 固定体積への切替後の確認 |
| `prod_rt` | NVT、300 K | 最初20 ns | 20,000,000 | 室温の構造・運動の解析 |

基本経路の合計は**67.4 ns＋最小化**。平衡化不足なら延長する。0.2 nsの加熱時間も最初の設計値であり、不安定なら加熱を緩やかにする。NPTの100 barは参照に合わせた準備操作で、室温物性を100 barで評価する意味ではない。

`eq_rt`の最終箱が採用解析区間の平均体積から1%以上外れている場合、固定体積Productionの代表性を再検討し、NPTを延長するなどして判断する。NVTでは体積が一定なので、**`prod_rt`の密度が横ばいであることを再平衡化の証拠には使わない**。NVTの平均圧力が1 barから有意に外れる場合も、体積と事前NPTを点検する。

参照との差分は、PMMA/OTfへの分子変更、力場・電荷の新規作成、PMEを使う最小化、初回の緩やかな加熱、段間状態の継続、追加冷却・室温段、室温ProductionのNVT化である。参照の軌跡そのものの再現とは呼ばない。

### 9.2 共通MD設定の設計値

次は本教材の初期設定である。MDP生成時には単位、適用するstage、実機で解決された値を確認し、重複キーを禁止する。

| 分類 | 基本値 |
|---|---|
| 積分 | `integrator=md`、`dt=0.001 ps`、`pbc=xyz` |
| 非結合 | `cutoff-scheme=Verlet`、`nstlist=20`、`verlet-buffer-tolerance=0.005` |
| 静電 | `coulombtype=PME`、`rcoulomb=1.0 nm`、`ewald-rtol=1e-5` |
| PME格子 | `pme-order=4`、`fourierspacing=0.12 nm`。実際の格子寸法も保存 |
| LJ | `vdwtype=Cut-off`、`rvdw=1.0 nm`、`vdw-modifier=Potential-shift`、`DispCorr=EnerPres` |
| 拘束 | `constraints=h-bonds`、`constraint-algorithm=lincs`、`lincs-order=4`、`lincs-iter=1` |
| 温度 | `tcoupl=v-rescale`、`tc-grps=System`、`tau-t=1.0 ps` |
| NPT | `pcoupl=c-rescale`、`pcoupltype=isotropic`、`tau-p=5.0 ps`、`compressibility=4.5e-5 bar^-1` |
| NVT | `pcoupl=no`。圧力制御パラメータをNPTから誤継承しない |
| COM | `comm-mode=Linear`、`comm-grps=System`、`nstcomm=100` |
| エネルギー | `nstcalcenergy=100`、`nstenergy=1000`、`nstlog=1000`（記録は1 psごと） |
| 座標 | `nstxout=0`、`nstvout=0`、`nstfout=0`、`nstxout-compressed=1000`、全原子 |
| XTC精度 | `compressed-x-precision=10000`（0.0001 nm）。短距離・低移動量の分解能も検査 |
| 初速度 | `heat`だけ`gen-vel=yes`、`gen-temp=300 K`、`gen-seed=101` |
| 段間継続 | それ以降`gen-vel=no`、`continuation=yes`、前段cptを使用 |
| stochastic seed | 基本replicaは`ld-seed=201`。実際に採用されたseedを保存 |
| 冷却以外 | `annealing=no`を明示。`heat`のみ加熱、`cool`のみ冷却 |

`heat`はMD最初の段なので`continuation=no`とし、最小化後にh-bondsを適用する。`annealing=single`、2点の時刻0/200 ps、温度300/1400 Kで加熱する。この段だけ速度を生成し、後段では前段の速度を継承する。

PMMAの圧縮率として水のような値が正しいと主張しない。ここではバロスタットの初期制御値であり、平均物性のfitには使わない。体積振動・不安定があればτPと併せて検査する。種別ごとのCOM除去や、小さいイオングループごとの過剰な温度制御は行わない。

温度・圧力の目標値は各stageから取得する。最小箱長は、全NPT軌跡について少なくとも最大非結合カットオフの2倍を超えることを検査し、さらに実効neighbor-list半径の制約も満たす。`rlist`等の実効値はgrompp出力から記録する。

上表は出発点であり、高温の時間刻み・拘束・cutoffの妥当性は短いpreflightと必要な比較計算で検査する。MDPキーの意味は[GROMACS 2026.3公式パラメータ説明](https://manual.gromacs.org/documentation/2026.3/user-guide/mdp-options.html)を参照する。

### 9.3 時間・checkpointの接続を曖昧にしない

MD開始を0 psとする累積時間をmanifestで管理する。`tinit=0`、`init-step=前段の最終global step`、`nsteps=今回追加するstep数`を生成し、前段checkpointのstep・timeと照合する。EMのstepはMDの累積時間に加えない。dtを変えた場合はこの単純なstep換算を使わず、接続時刻を再設計する。

延長がない基本経路の累積時刻は、`heat:0–200`、`eq_1:200–700`、`eq_2:700–1200`、`eq_3:1200–6200`、`prod:6200–26200`、`cool:26200–37200`、`eq_rt:37200–47200`、`nvt_rt:47200–47400`、`prod_rt:47400–67400 ps`である。

段を変えるときは`grompp -t previous.cpt`で座標・速度を渡し、必要な制御状態も引き継ぐため`mdrun -cpi previous.cpt`を使う。アンサンブルや制御条件の変更点は記録し、変更後の段で安定を再確認する。[grompp公式説明](https://manual.gromacs.org/documentation/2026.3/onlinehelp/gmx-grompp.html)に沿って実装し、**gromppの`-t`だけで完全な状態継続と断定しない**。

例として、`eq_2`の作業ディレクトリで、前段の実ファイルをmanifestから解決して実行する。

```bash
gmx grompp -f eq_2.mdp -c ../eq_1/eq_1.gro -t ../eq_1/eq_1.cpt -p topol.top -o eq_2.tpr -po eq_2_effective.mdp -pp eq_2_processed.top
gmx mdrun -deffnm eq_2 -cpi ../eq_1/eq_1.cpt -noappend -cpt 5
```

`-noappend`では`.partXXXX`付き出力になる場合がある。stage表の`eq_2.gro`等は論理名であり、**次段へ渡す実ファイルはmanifestから取得**する。実装ではディレクトリの存在しない例示パスをそのまま呼ばない。partファイルの出現・最終出力名・global stepを、短い2段接続テストで確認してから長時間計算する。

同一stageを中断再開するときは、新しい初速度を作らず、同じTPR・同じ出力名で`mdrun -deffnm <stage> -cpi <checkpoint> -append`を使う。出力checksumが一致しないときは原因を調べ、正本を復元するか、別partへ継続して解析で結合する。出力を消して警告を回避しない。[再開と出力追記の公式説明](https://manual.gromacs.org/documentation/2026.3/user-guide/managing-simulations.html)を参照する。

**接続テストの必須検査：** 前段cptの時刻→新TPRの開始時刻→新ログの開始/終了stepが一致すること、要求した追加時間だけ進むこと、速度再生成がないこと、変更したT/P設定が有効なこと、冷却目標が開始時点で1400 Kになっていること。

## 10. Phase 4：最小化と短時間preflight

### 10.1 最小化

`opt.mdp`は`integrator=steep`、`nsteps=50000`、`emtol=100 kJ mol^-1 nm^-1`、`emstep=0.01 nm`、温度・圧力制御なし、速度生成なしとする。PME、1.0 nmの非結合条件を使い、本計算と大きく異なる静電モデルで最小化しない。EMの結合拘束は基本`constraints=none`とし、MD開始時にh-bondsを適用する。

grompp後に0-step評価を行い、NaN/Inf・異常な近接・過大な力の原子を特定してからEMを実行する。初期のEpotが正か負かだけで判定しない。

**合格基準：** 原則Fmax≤100、終了理由と最大力原子を保存。machine precisionで止まり目標未達なら「目標達成」と書かない。重なり、力場、力が集中した箇所、step幅を点検し、改善後も未達ならFmaxと短時間安定性の根拠を示して明示的に条件付き判定する。参照Fmax≈231を無条件の許容基準にはしない。

### 10.2 MD preflight

本番とは別runで、同じ力場による10–50 ps程度の短いMDと、短い加熱・NPT・2段接続・中断再開を検査する。必要なら初期1–10 psを0.5 fsで検査し、1 fsとの差を調べる。小さい箱で行う動作テストは構成と密度を記録し、基本系の物性結果に混ぜない。

**合格基準：** Fatal、LINCS警告、NaN、温度暴走、極端な結合角の崩壊、異常な箱変化がない。熱浴のあるNPTで総エネルギー一定を要求しない。数値積分が疑わしい場合は、小時間刻み比較や短いNVE診断を用いる。保存量のドリフトだけを普遍的な単一閾値で合否判定しない。

**成果物：** `04_minimization_and_preflight.ipynb`、`opt.*`、短時間ログ、終了理由、構造の前後表示、preflight判定。

## 11. Phase 5：高温平衡化と高温Production

`heat → eq_1 → eq_2 → eq_3`を順に実行する。各段の温度・圧力・密度・体積・Potential・全stageの時間軸を可視化する。各段の時間が予定に達しただけでは平衡化合格にしない。

### 11.1 平衡化の判定方法

以下は教材の運用目安であり、普遍的な物理法則として扱わない。

- 完了予定stepとログの最終step、読み込めるedr/xtc/cpt、原子数、時刻の単調性を検査する。
- 定温区間で平均温度が目標の概ね1%以内であることを確認する。持続的な偏りなら原因を調べる。
- 圧力は瞬間値ではなくブロック平均・不確かさで評価する。1 barに対する相対誤差1%等を要求しない。
- まず緩和部分を除いた候補区間を選び、自己相関とブロック長依存を調べる。原則4ブロック以上、各ブロックを主要観測量の相関時間より十分長くする。
- 候補区間の前半/後半密度の差が概ね1%以下か、差の不確かさで整合するかを調べる。**不確かさが大きすぎるために差を検出できない場合は合格にしない**。
- Potentialもブロック平均と時間傾向を確認する。Rg、主鎖二面角、イオン接触数に一方向の変化が残っていないかを併せて見る。
- 密度が横ばいでも鎖・イオンの配置が十分サンプルされたとは限らない。独立replicaや追加時間で検査する。

条件を満たさない場合は、同条件で5 ns単位などの延長を提案し、延長理由・予算・判定区間を記録する。高圧の短い`eq_1/2`に最終構造統計の収束までは要求しないが、数値不安定な状態を次段へ渡さない。

### 11.2 高温Production

`eq_3`の判定後に`prod`を20 ns実行し、後述の標準解析を一度通す。高温での密度・移動・イオン会合を室温との比較用に保存する。

`prod`はNPTなので、長時間MSDのunwrapと箱変動の扱いを検査する。高温の拡散係数を定量比較する場合は、`prod`の採用状態から固定箱NVTの短い再平衡化とProductionを別枝として追加する。**NPTでの単純unwrapに由来する見かけの拡散を、そのまま室温NVTと比較しない**。この追加枝の時間・費用は計画に反映する。

**成果物：** `05_high_temperature_md.ipynb`、`eq_1/2/3`と`prod`の入力・出力、性能測定、平衡化判定、解析区間。

## 12. Phase 6：冷却・室温NPT再平衡化・室温Production

### 12.1 冷却条件

基本冷却速度を100 K/nsとする。これは分子シミュレーションで実行可能な学習用速度であり、実験の冷却速度の再現ではない。

`cool`はNPT / 1 bar / 11 ns、単一System温度群、`annealing=single`とする。相対時刻から次の点列を生成する。

| 冷却開始からの時刻 [ps] | 0 | 2000 | 4000 | 6000 | 8000 | 9000 | 9500 | 10000 | 10500 | 11000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 目標温度 [K] | 1400 | 1200 | 1000 | 800 | 600 | 500 | 450 | 400 | 350 | 300 |

実際の`annealing-time`は、manifestの冷却開始時刻を上表に足した**絶対シミュレーション時刻**で書く。基本経路での例：

```text
; cool.mdp の段固有部分。共通設定と合わせて生成する
tinit             = 0
init-step         = 26200000
nsteps            = 11000000
ref-t             = 1400
ref-p             = 1.0
gen-vel           = no
continuation      = yes
annealing         = single
annealing-npoints = 10
annealing-time    = 26200 28200 30200 32200 34200 35200 35700 36200 36700 37200
annealing-temp    = 1400 1200 1000 800 600 500 450 400 350 300
```

上流を延長したら、開始stepと時刻列を**自動再生成**する。時刻0開始の点列を累積時刻のrunへそのまま渡す事故を防ぐ。preflightでは目標温度列と実温度の追従を必ず図示する。冷却中は非平衡なので、この全区間をまとめたMSDから単一のDをfitしない。

### 12.2 室温再平衡化

室温は300 Kと定義する。`eq_rt`は`annealing=no`、300 K / 1 bar NPTでまず10 ns。体積・密度・Potential・Rg・二面角・接触数を検査し、Phase 5の方法で安定性を判定する。

PMMA系が室温でガラス状になり、密度や局所構造は安定しても大域的な平衡に達しない可能性がある。その場合は「この冷却履歴に対応する準安定ガラスモデル」「観測時間内の局所安定」として報告する。時間を延ばせば必ず平衡になると約束しない。

`eq_rt`後、最終箱を固定して`nvt_rt`を0.2 ns、続いて`prod_rt`をまず20 ns行う。室温構造・運動の解析には後者の定温区間を使う。前半/後半でPotential・接触数・局所運動が変わるなら、エージングを疑い、独立に報告する。

**成果物：** `06_cooling_and_room_temperature.ipynb`、冷却のT(t)/ρ(t)/V(T)、`eq_rt.*` / `nvt_rt.*` / `prod_rt.*`、室温状態の判定と延長判断。

**合格基準：** 300 Kに到達、annealingが後段で無効、室温のT・体積・構造の判定根拠あり。準安定状態は準安定として区別されている。拡散が小さいことだけをMD失敗としない。

## 13. Phase 7：標準解析と結果の検証

### 13.1 解析データと単位

各解析は元tpr/xtc/edr、使用区間、間引き、対象atom/molecule、前処理、コード版を`analysis_config.yaml`へ保存する。Notebookを上から実行すればCSVと図が生成されること。手入力で結論を埋めない。

GROMACSはnm / ps / kJ mol⁻¹ / bar / K、MDAnalysisは通常Å / psである。内部解析単位はnm / psに統一し、明示的に変換する。`1 Å²=0.01 nm²`、`1 nm²/ps=1e-6 m²/s`、密度は`kg/m³÷1000=g/cm³`。

解析用に2系統の座標を作る。

- **構造・RDF用：** 分子をつなぎ、各frameのboxに対応した最小イメージ距離を使う。
- **MSD用：** 各分子をwholeにした後、分子COMの連続性を保ってunwrapする。表示用のcenteringや`-pbc mol`だけでMSD前処理が完了したと思わない。回転fitや鎖ごとの重心移動除去で並進を消さない。

既知の箱横断と既知の並進を持つ小さな検証軌跡で、PBC前処理・COM・単位を検査する。NPTでは箱変動に対応したunwrap方法を別途検証する。保存間隔の間に半箱を超える移動があり得るなら、座標の跳びからの再構成は不定になるため保存間隔を短くする。

### 13.2 密度・エネルギー

密度の独立検算式は次を使う（Mは箱中の全原子質量の合計、単位amu）。

```text
ρ [g/cm³] = 1.66053906660e-3 × M [amu] / V [nm³]
V = |det(box_vectors)|
```

edrはpyedrと`gmx energy`で同一区間の値を比較する。少なくともTemperature、Pressure、Potential、Kinetic、Total Energy、Volume、Density、箱長を出力する。存在する保存量も診断用に記録する。

密度・Potential・圧力は、ブロック平均・不確かさ・前半/後半差・採用区間を表にする。Potentialの絶対値を、化学種数や力場の違う参照系と直接比較しない。必要なら1原子・1質量あたり等の規格化も併記する。NVTの固定密度は状態設定であり、独立に予測した室温密度は`eq_rt`のNPTから報告する。

### 13.3 MSDと拡散係数

各分子COMについて、複数の時間原点を用いて計算する。

```text
MSD(τ) = average over molecules i and valid origins t0 of
         |R_i(t0+τ) - R_i(t0)|²
D = (1/6) × slope of MSD(τ)       # 等方3次元・Fick拡散の場合
α(τ) = d log(MSD) / d log(τ)
```

TPS、OTf、PMMA鎖COMを分ける。PMMAのセグメントMSDは鎖COMの並進と別に評価する。fitに切片を許すかも記録する。

必須の図は、MSDの線形軸、log-log、局所傾きα、fit区間、時間原点数または各lagの有効サンプル数である。最長lagはまず解析区間長の1/4以下とし、相関した時間原点を独立サンプルとして数えない。

**D採用条件：** 数値・PBC検査に合格し、後期にα≈1の区間があり、fit範囲を変えてもDが大きく変わらず、ブロック/replica間で評価可能なこと。教材の診断目安はαが概ね0.9–1.1、fit範囲変更によるD差が概ね20%以内とする。ただしαのノイズと不確かさも確認し、この数値だけで自動合格にしない。

ケージ内振動、MSD plateau、subdiffusion、観測時間不足なら、`D: not_resolved`とする。自動fitが返す小さな傾きを採用せず、観測時間とMSD、plateau、α、移動量分布を報告する。`D=0`や厳密な上限を根拠なく宣言しない。

`gmx msd`による分子COM計算とも照合する。分子COMには`-mol`を用い、対象は分子全体が含まれるgroupとする。fit区間は必ず明示し、既定の10–90%を根拠なく採用しない。[GROMACS 2026.3のMSD説明](https://manual.gromacs.org/documentation/2026.3/onlinehelp/gmx-msd.html)も、線形区間外のDと誤差の解釈に注意を求めている。

有限サイズ補正は自動適用しない。適用に必要な粘度・流体の仮定などが未検証なら、その制限を明記した箱サイズ依存の値として報告する。

### 13.4 RDF・配位数

必須RDFを次とする。解析サイトは`site_map.json`から取得する。

| 対象 | わかること | 解釈の注意 |
|---|---|---|
| TPS中心S — OTf S | イオン間の代表距離。参照との解析形式の対応 | 別アニオンなのでピーク位置・高さの一致を要求しない |
| TPS中心S — OTf O | カチオンとスルホン酸基の接近 | O配位数とアニオン個数は異なる |
| TPS中心S — TPS中心S | カチオン同士の配置相関 | 同じ分子自身を除く |
| OTf S — OTf S | アニオン同士の配置相関 | 同符号のピークだけで結合やクラスターを断定しない |
| TPS中心S — PMMAカルボニルO | ポリマーとの局所接触 | 特定の結合生成を意味しない |
| OTf S — PMMAカルボニルO | アニオン周囲のポリマー配置 | 元素Sの部分電荷とアニオン全体の電荷を混同しない |

bin幅はまず0.01 nm、`rmax`は解析区間中の最小箱長の半分より小さくする。規格化、同一種self除外、分子内除外、有限粒子数の扱い、frameごとの体積の扱いを明記する。自作RDFと`gmx rdf`の同一種規格化が異なる場合は、その差を説明してから比較する。[公式RDF説明](https://manual.gromacs.org/documentation/2026.3/onlinehelp/gmx-rdf.html)ではexclusionを指定しても規格化の粒子数は自動変更されない。

配位数は`4πρ_B ∫g_AB(r)r²dr`と、frameごとの直接近接数の両方から検算する。NPTでは体積変動の平均方法が一致するか確認する。

接触距離は各RDFの第一極小から決める。極小が曖昧なら、固定した候補距離と±0.05 nm等での感度を示す。参照の6.5 Åを別イオンへそのまま転用しない。高温/室温の比較では共通cutoffでの比較を基本とし、条件別第一極小を使う結果は別表示する。

**成果物：** `07_standard_analysis.ipynb`、時系列・MSD・RDF・配位数のCSVと図、fit記録、D採用/未決定の判定表。

## 14. Phase 8：レジスト関連の分子挙動を解析する

### 14.1 必須解析

| 問い | 計算する量 | 必須の出力と検証 |
|---|---|---|
| 冷却でイオン会合はどう変わるか | TPS–OTf接触率、free ion率、距離分布 | 高温/室温の比較、cutoff感度、誤差 |
| イオンクラスターはあるか | 接触グラフの連結成分、サイズ・組成・正味電荷 | 単量体/1対/多体の割合、数平均と粒子重みの区別 |
| 接触がどれくらい持続するか | イオン対の連続/間欠接触相関 | 生存曲線、イベント数、観測窓、打切りの扱い |
| イオンがポリマーに接触するか | カルボニルOへの接触数、近接鎖数、滞在相関 | Oの個数と繰り返し単位/鎖の個数を区別 |
| 鎖の大きさが変わるか | 分子whole後のRg、末端間距離、分布 | 鎖ごとの分布・時間変化、PBC破綻の除外 |
| 局所運動が遅くなるか | 主鎖二面角分布/遷移、セグメントMSD、配向相関 | 原子定義、サンプル数、高温と室温の比較 |
| 移動するイオンと動かないイオンがあるか | 時間窓別のCOM変位分布、接触状態別MSD | window依存、状態選択バイアス、少数事象の区別 |

接触グラフの辺は、**1つのTPS中心Sと1つのOTfの3つのOとの最短距離**が採用cutoff未満であること、と定義する。Oが2個近くても同じイオン対の辺は1本と数える。S–S RDFの距離と、このS–O接触cutoffを混同しない。ポリマーへの接触も、近接原子数・単位数・鎖数を分ける。

連続相関は接触が一度も途切れない確率、間欠相関は途中の解離・再会合を許す定義とし、式・時間分解能・短い途切れの許容の有無を保存する。観測終了時まで残る接触は打切りとして扱う。数件以下の解離しかない場合、信頼できる単一寿命を出さず、イベント数・生存曲線・観測時間で報告する。

高温から室温で同じ接触が増えたことだけでは、必ずしも熱力学的な結合強化や平衡定数を示さない。冷却中の凍結、初期状態、有限サイズ、電荷モデルへの依存を検討する。実験値を比較する場合は、分子量・配合・温度・測定条件の合う一次資料を別途確認する。

### 14.2 発展解析（基本到達点と分離）

以下は追加目的に応じた別計画とする。

- 冷却速度50/100/200 K/ns等の比較。密度・接触・移動の熱履歴依存を調べる。
- 温度ごとの保持計算を追加し、V(T)の傾き変化を検討する。単一の急冷曲線から実験Tgを断定しない。
- probe半径・原子半径・格子幅を固定した自由体積解析。単に低密度部分を自由体積と呼ばない。
- 複数の立体配列、鎖長、PAG配合、箱サイズの比較。
- 検証済み別電荷モデルとの比較。電荷を任意に縮めて拡散を合わせない。
- 反応・酸・表面を含む研究へ拡張する場合は、新しいモデルの仕様を別途作る。

**成果物：** `08_resist_molecular_behaviors.ipynb`、接触定義、クラスター・寿命・鎖/局所運動のCSVと図、各問いへの短い結論と制限。

## 15. Phase 9：再現性・統計・最終報告

### 15.1 replicaと完了の水準

まず1本で一連の操作を習熟し、その後**独立初期配置から基本経路を3本**実行して、結論の配置依存を調べる。同じcheckpointから速度だけ変えた枝は、独立な初期構造replicaと区別する。

基本の立体配列と力場は固定し、例として次のseedをconfigに明記する。

| replica | PMMA/TPS/OTf挿入seed | 初速度seed | thermostat等のseed |
|---|---|---:|---:|
| r01 | 1 / 2 / 3 | 101 | 201 |
| r02 | 11 / 12 / 13 | 111 | 211 |
| r03 | 21 / 22 / 23 | 121 | 221 |

frame数を独立サンプル数としない。時系列内のブロック解析と、replica間のばらつきを分けて報告する。3本は十分性を自動保証せず、結論が揃わない場合は不確かさまたは追加計算の必要性を示す。

完了表示を次の3段階に分ける。

- **実装・実行習熟完了：** r01の全経路と必須解析を実行し、各ゲートの根拠と確認問題への回答がある。
- **基本比較完了：** r01–r03を集計し、高温/室温の結論とばらつきを報告した。
- **研究用妥当性検討：** 必要な力場検証、時間・サイズ・冷却履歴依存、実験/文献比較を追加している。基本比較完了と同義にしない。

### 15.2 最終成果物

`09_validation_and_report.ipynb`から、最終Markdown報告と図・CSVを生成する。PDF化は任意とする。報告には以下を含める。

1. 分子構造、PMMA末端・立体配列・重合度、分子数、配合、力場・電荷モデル。
2. 実環境、入力hash、全stageのT/P/時間・実行時間・終了判定、再開履歴。
3. NPTでの高温・室温密度と不確かさ、体積・温度・エネルギーの診断。
4. 冷却履歴と室温での安定性。平衡/準安定/エージングの判定。
5. TPS・OTf・PMMAのMSD、Dまたは`not_resolved`、fit区間・単位・不確かさ。
6. RDF、配位数、イオン対・クラスター・接触寿命・ポリマー運動。
7. 独立replica間のばらつき、支持された結論、未解決事項、モデルの限界。
8. 元入力から再実行するコマンドと、解析だけ再実行するコマンド。

## 16. PBS実行・時間予算・復旧の要件

### 16.1 実行スクリプト

現在の成功済みPBSジョブを基に作り、既存リポジトリのSlurmテンプレートをそのまま使わない。`ncpus=16`という要求だけではMPI16 rankとは決まらないため、launcherとGROMACSログの実rank/threadを照合する。

- 原則1stage/1job。依存投入は前stageの合格判定後に行う。
- `set -euo pipefail`等を使い、gromppや解析前検査の失敗を次段へ伝播させる。
- PATHと環境有効化をジョブ内で明示し、対話shellで見えたgmxがPBSでも同じか確認する。
- OpenMP/MPI/GPU設定は成功済み設定を出発点とし、実測なしで最適と呼ばない。
- walltime内に収まるよう`-maxh`を設定し、終了余裕を残す。checkpointは例として`-cpt 5`で5分ごと。
- `-reprod`は比較実験で使用可。ただしバージョン・rank等が異なる軌跡の完全一致を保証しない。
- ジョブID、投入・開始・終了、return code、実資源、ns/day、stdout/stderrをmanifestに保存する。
- `grompp -maxwarn`で未知の警告を一律回避しない。既存エタノール教材の`comm-mode=none`例外はこの系へ持ち込まない。

### 16.2 所要時間と容量

参照の約23.7 ns/dayは先輩の別組成・別バージョンでの性能であり、予定値にしない。基本系の0.1–0.5 ns程度で実測し、次の式で予算を出す。

```text
walltime [day] = planned trajectory [ns] / measured performance [ns/day]
```

仮に23.7 ns/dayなら基本経路67.4 nsは約2.84日/replica、3本で約8.53計算日相当になる。これは計算例であり、QM/力場準備、queue待ち、延長、解析・追加枝は含まない。30時間walltimeの1jobに全経路を詰めない。

1 ps間隔で保存する67.4 nsは約67,400 frame/replicaとなる。参照より高頻度・高精度にするため、容量は短時間出力から外挿する。edr/log/cpt/一時TRR/解析用軌跡を含め、見積りの2倍程度を目安に空き容量を確保する。XTCサイズを原子数だけから確定しない。容量を減らすなら、最初に必要な時間分解能を検討する。

Notebookの「全セル実行」で長時間MDや既存runの上書きを始めない。準備・解析とジョブ投入を分け、`run_id`の新規性、stage指定、完了済みrunへの書込み禁止を実装する。投入済みjob IDがある場合は再投入前に状態を照会する。

### 16.3 失敗時の分岐

| 症状 | 調べるもの | 対応 |
|---|---|---|
| 原子数・原子名が不一致 | atom-map、変換、分子順序、H | 入力を修正。警告の抑制で進めない |
| 総電荷が非ゼロ | 形式電荷、mol2→itp、丸め、分子個数 | 元データから修正。PME背景電荷で妥協しない |
| 力場項が未定義 | 原子型、frcmod、S/主鎖二面角 | 出典・検証を追加してから再開 |
| 挿入数が不足 | 配座、箱、半径、実挿入数 | 条件を見直し、configと密度を更新 |
| EM未収束 | 最大力原子、重なり、パラメータ、終了理由 | 原因別に修正または条件付き判定を記録 |
| LINCS/温度暴走 | 近接、分子内変形、dt、加熱、拘束 | 直前の健全な状態から短時間検査。警告閾値を緩めない |
| NPT体積が暴れる | 初期圧力、密度、τP、compressibility | 制御・構造を検査。prodへ進めない |
| stage接続で時刻が不整合 | init-step、nsteps、tinit、cpt、part | 2段接続テストへ戻り、時間設定を再生成 |
| 室温でMSD plateau | α、時間窓、粒子分布、PBC、XTC精度 | 未観測/非Fickとして記録。Dを強制fitしない |
| 接触の解離が少ない | event数、観測窓、打切り | 寿命未決定と生存曲線を報告 |
| walltime終了 | job状態、最後のcpt、出力checksum | 同じstageを継続。自動で最初から再投入しない |
| ログだけ終了・軌跡破損 | gmx check、予定step、容量、I/O | 未完了として復旧し、完走扱いしない |

## 17. 習熟Notebookの必須構成

全Notebookを、日本語で次の順に構成する。処理の中核はscripts/moduleへ置き、Notebookはその処理と結果の理解に集中する。

1. **今回の問いと前提：** 何を知りたいか、どの入力から始めるか。
2. **物理的な説明と式：** 例えば密度なら「質量÷体積」、MSDなら「一定時間後の移動距離の二乗平均」。変数と単位を説明する。
3. **Python/コマンド：** 式のどの部分を計算するかを対応付ける。乱数・入出力先を明示する。
4. **可視化：** 軸・単位・区間・対象分子・条件を図に書く。
5. **観察：** 実際の出力から短い説明を作る。期待した傾向と違っても結果を改変しない。
6. **合格判定と限界：** 完走・安定・統計・モデル妥当性を区別する。
7. **確認問題と次段の入力：** 自分で答えられることを確認する。

代表的な確認問題：

- なぜPMMA19を選び、末端Hが分子式にどう入るのか。
- なぜF1のitpをOTfへ流用できないか。
- 100鎖と83イオン対から、総電荷・原子数・塩の質量分率を計算できるか。
- `[ defaults ]`と`[ atomtypes ]`をどこへ置くか。1–4の倍率は何に作用するか。
- 最小化のstepとMDの時間はどう違うか。Fmax目標未達をどう判断するか。
- 速度生成、thermostat、checkpoint再開はそれぞれ何をするか。
- 1400 Kと300 Kで観測できること・解釈できないことは何か。
- NPT密度が横ばいでもポリマーが平衡とは限らないのはなぜか。
- NVTで密度が横ばいなのは何を意味するか。
- MSDが曲線やplateauなら、なぜ直線fitの数値をDに採用できないか。
- RDFのO配位数とイオン対数はどう違うか。
- 「OTfが遅い」と「酸の拡散が遅い」は、なぜ同じ結論にならないか。

## 18. Codexの実装契約・検証・進捗台帳

### 18.1 Phase単位の作業

Codexは、対象Phaseの実ファイルと前Phaseの判定を読み、今回の範囲を短く説明してから実装する。**Phase内の低リスクな修正と検証はまとめて完了**し、コマンド1つごとの確認を挟まない。

利用者が「Phase Nまで」と指定したら、その範囲を越えない。未承認の環境変更、モデル/配合/力場の変更、予算を大きく増やす延長、長時間ジョブの投入は、具体的な差分と予算を用意して判断を求める。あらかじめ一連の実行が承認されている場合は、その範囲で繰り返し承認を求めない。

各Phase終了時に、変更ファイル、実施した検査、実行結果、未検証点、合格/条件付き/不合格、次に使う入力を`PROGRESS.md`へ記録する。

| 状態 | 意味 |
|---|---|
| planned | 仕様のみ |
| implemented | 実装済み、対象環境で未実行 |
| running | job ID等の実行証拠あり |
| executed | 出力生成済み、判定未完了 |
| passed | このPhaseの必要な検証を合格 |
| conditional | 目標未達や限定条件と、その根拠を明記 |
| failed | 問題を検出、次段への入力として未採用 |

manifestには最低限、run ID、replica、親stage、入力/出力hash、git commit＋未commit差分の記録、config、実コマンド、tool実体/version、seed、開始/終了時刻、MDのstep/time、PBS job ID、実出力パス、判定を持たせる。予定値と実測値を同じ欄へ上書きしない。

### 18.2 必要な検証

実装時には、科学的な誤りを防ぐ次の検査を行う。長時間MDをユニットテストのたびに実行しない。

- PMMA結合グラフ、末端・式・質量・立体配置、TPS/OTfの電荷と原子対応。
- 力場項・型衝突・1–4・電荷合計・トポロジーと座標の順序。
- `dt×nsteps`、累積時刻、冷却点列、速度生成回数、NVT/NPT設定、MDP重複キー。
- 既知の箱横断・分子COM・距離・密度・RDF・既知変位のMSDに対する解析検査。
- 実GROMACSでのgrompp、0-step、短い2stage接続とcheckpoint中断再開。
- Notebookのnbformat、コード構文、指定カーネルでの該当Phase実行、保存出力・図の確認。
- 共有moduleを変更した場合は影響する既存テストを実行し、既存エタノール系の動作境界を維持する。

`gmx check`で軌跡の原子数・時間・読込を確認し、最終stepはログ/checkpointと照合する。出力間隔によって最終stepの座標がXTCにない場合を、正常な最終gro/cptと区別して判断する。複数partの結合時は重複frame・時刻の飛び・原子順を検査し、元partを保存する。

### 18.3 開始時にCodexへ渡す指示

```text
MD習熟/PMMA_TPS_OTf_MD_tutorial_spec.md を正本として、Phase 0を実装してください。
参照Notebookと既存コードを確認し、再現成功環境の実体を確認できる形で環境検査、
config、進捗台帳、00_environment_and_plan.ipynbを作ってください。
基本案はPMMA19(H終端)100鎖、TPS+83、OTf-83、GAFF2+AM1-BCCです。
仕様上の設計値と実測値を区別してください。
今回は環境の新規導入・長時間ジョブ投入・Phase 1以降の実装は行わず、
検証結果と不足情報を示してPhase 0の範囲で止めてください。
```

以後は次を用いる。実行を依頼する場合には、許可するstageと時間予算も指定する。

```text
仕様書とPROGRESS.mdを読んで、次のPhase Nを実装・検証してください。
前段の合格証拠と実入力を確認し、Notebookは物理の説明→式→Python→図→解釈の順にしてください。
今の承認範囲はPhase Nまでです。実施内容、実行/未実行、合格基準の判定を記録してください。
モデル変更や長時間実行が必要なら、具体的な変更案・費用と判断点を用意してください。
```

## 19. 受入チェックリスト

- [ ] 実行環境が利用者の再現成功環境と対応し、実体・版・PBS設定を記録した。
- [ ] 3種を一から定義し、PMMA19の末端・立体配列・原子数・式・分子量を検査した。
- [ ] 電荷モデル・力場・変換過程と未解決パラメータの有無を確認した。
- [ ] 100/83/83、32,186原子、中性、塩約15.24 wt%を実topologyから検算した。
- [ ] 挿入・最小化・preflight・2stage接続・checkpoint再開を検証した。
- [ ] `eq_1 → eq_2 → eq_3 → prod`を実行し、数値安定・構造安定を判定した。
- [ ] 冷却の目標温度と実温度、時刻接続を確認した。
- [ ] 300 KでNPT再平衡化と固定体積Productionを行い、状態の限界を報告した。
- [ ] 密度・エネルギー・MSD・D評価可否・RDF・配位数を根拠付きで示した。
- [ ] イオン対・クラスター・接触寿命・ポリマーとの接触・鎖/局所運動を解析した。
- [ ] 1本の習熟完了と3本の基本比較完了を区別し、必要な段階まで実施した。
- [ ] 全Notebookに説明・式・コード・図・確認問題・合格判定がある。
- [ ] 再実行手順、入力hash、ログ、checkpoint、解析条件、最終報告を保存した。
- [ ] 高温処理・準安定ガラス・固定電荷・非反応性・短鎖・有限サイズの制限を書いた。

## 20. 参照資料

参照確認日：2026-10-27。本書で外部仕様の根拠とした一次資料を示す。リンクの最新版へ自動追従せず、実機2026.3の`-h`と照合する。

- [このリポジトリの参照Notebook](MD_output_guide.ipynb)：構成、元系の分子量、T/P/時間、保存出力、PBS実行条件。
- [このPCの環境記録](../docs/current_pc_setup.md)：Mac側環境の既存記録。HPCの現状とは区別する。
- [GAFFの公式解説](https://ambermd.org/antechamber/gaff.html)：力場と電荷モデルの背景。GAFF2の実データは採用AmberTools内で固定する。
- [Antechamberの公式解説](https://ambermd.org/antechamber/ac.html)：構造・型・電荷・パラメータ生成。古い例示もあるため現行CLIは実機で検査する。
- [GROMACS 2026.3：MDP options](https://manual.gromacs.org/documentation/2026.3/user-guide/mdp-options.html)：MDP・annealing・時間・拘束・制御。
- [GROMACS 2026.3：insert-molecules](https://manual.gromacs.org/documentation/2026.3/onlinehelp/gmx-insert-molecules.html)：挿入方法、半径、seed。
- [GROMACS 2026.3：grompp](https://manual.gromacs.org/documentation/2026.3/onlinehelp/gmx-grompp.html)：入力処理、checkpoint入力、実効MDP。
- [GROMACS 2026.3：Managing long simulations](https://manual.gromacs.org/documentation/2026.3/user-guide/managing-simulations.html)：再開、追記、part出力、延長。
- [GROMACS 2026.3：msd](https://manual.gromacs.org/documentation/2026.3/onlinehelp/gmx-msd.html)：分子COM・fit・線形性・時間原点。
- [GROMACS 2026.3：rdf](https://manual.gromacs.org/documentation/2026.3/onlinehelp/gmx-rdf.html)：RDF規格化、exclusion、配位数。
