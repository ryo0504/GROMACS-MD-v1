# PMMA / TPS / OTf：研究室サーバーへ持っていく入力セット

準備日：2026-10-28。実行予定：2026-10-29（JST）。

このフォルダは、[v2仕様書の保存コピー](references/PMMA_TPS_OTf_MD_tutorial_spec_v2.md)に基づく**入力準備セット**です。分子構造の作成と軽い形式検証は行いましたが、実分子のGaussian・ESP・RESP・力場変換・MDは実行していません。保存コピー内の相対パスは元repoの配置を指します。今回の実行条件は本README・decisions.md・configを優先してください。

最初に[明日の実行手順](docs/01_lab_server_steps.md)、仕組みを理解するときは[電荷からMDへつなぐ解説](docs/02_charges_to_md.md)を読んでください。

## 確定した前提

- MDは再現計算に使った研究室のLinux/PBS環境で行います。Macでは重い計算をしません。
- 分子はH終端PMMA19（100鎖）、TPS⁺（83個）、CF₃SO₃⁻（83個）。初期箱は10 nm角、中性、32,186原子です。
- 力場はGAFF2＋RESP。電荷scaleは1.0です。実計算前に電荷や結合定数を仮の値で埋めません。
- **Gaussian 16、RHF/6-31G*、MK格子ESP＋AmberTools RESP**を採用します。2026-10-28の利用者回答による決定です。
- v2のPySCF・独自Fibonacci格子から、Gaussian 16・MK格子へ変更した点は[decisions.md](decisions.md)へ記録しました。元のMDの1400 K、圧力、時間、段名は維持します。
- Gaussianに使えるCPU・メモリ・PBS queue・環境有効化方法は不明です。`config/server.example.json`は未設定で、**そのままではジョブを生成できません**。明日確認します。

## 何が入っているか

| 場所 | 今回作成したもの | 状態 |
|---|---|---|
| `structures/` | PMMA19/3/5、TPS、OTfのSDF・XYZ・初期GRO、原子対応JSON/CSV | 初期3D構造。QM最適化済みではない |
| `gaussian/<分子>_cXX/` | `opt.gjf.in`と`esp.gjf.in` | 資源値を設定すれば実行入力へ展開できる雛形 |
| `resp_templates/` | RESP拘束入力、内部/末端templateの変数と対応 | 原子対応・形式は検査済み。実ESPに合わせた重みは後で生成 |
| `config/` | 組成・計算条件・電荷プロトコル・サーバー設定例 | 正本の設定 |
| `mdp/` | opt、heat、eq_1/2/3、prod、cool、eq_rt、prod_rt、NVT追加枝 | v2のproduction時間に対応する入力 |
| `topology/topol.top.in` | 混合系topologyの構成と個数 | `.itp`完成前の雛形。実行可能な`.top`と区別 |
| `scripts/` | サーバー検査、Gaussian入力展開、ESP抽出、RESP、力場生成・変換検査、packing、stage準備 | 軽い検査済み。g16/実力場/MDのサーバー実行は未検証 |
| `tests/` | 構造、原子対応、MDP、RESP拘束・凍結、入力展開のテスト | 人工ESPはテスト専用 |
| `reports/` | 準備状態、検査結果、未実行の変換レビューゲート | 予定と実行済みを分けて記録 |

`@NPROC@`と`@MEM_GB@`のある`.gjf.in`を直接g16へ渡さないでください。`server.py render`で、確認した資源値を埋めた`.gjf`とPBSジョブを`runs/qm/`へ生成します。

## 分子とGaussianジョブの数

| 分子 | 原子数 | 電荷 | 用途 | 初期配座 |
|---|---:|---:|---|---:|
| PMMA19 | 287 | 0 | 実際に箱へ入れる鎖。全鎖QMは行わない | 1 |
| PMMA3 | 47 | 0 | 内部・両末端のRESP fitと評価 | 3 |
| PMMA5 | 77 | 0 | 同上。鎖長の影響を評価 | 3 |
| TPS | 34 | +1 | イオンRESP | 3 |
| OTF | 8 | −1 | イオンRESP。最初の確認ジョブ | 1 |

Gaussianジョブは合計10個です。1ジョブ内で構造最適化→ESP single point→ESP抽出を行います。PMMAのc01/c02をfit用、c03を評価用に分けています。最適化で配座が同じ極小へ収束することもあるため、最終配座が独立かどうかは実結果で検査します。

## 軽い初期構造処理について

RDKitで明示的H・結合・立体配置を持つ構造を作り、PMMAは初期座標調整としてMMFF94を最大300反復だけ使用しました。**PMMA19とPMMA5ではこの初期調整の収束条件に達していません**。これはGaussian最適化やGAFF2最小化の収束とは別です。PMMA3/5は明日Gaussianで最適化します。PMMA19は全鎖QMをせず、最終的にGAFF2で最小化します。

TPS/OTfはETKDGの初期配座です。ETKDG内部のUFF型判定で荷電Sに関する注意が出ましたが、UFFを本番力場として採用したわけではありません。結合長・原子間距離・指定立体配置は軽い検査を通しています。実QMでの構造とGAFF2のSパラメータは別途検証します。

## 最初に実行するコマンド

フォルダ全体を研究室サーバーへコピーして、適切なPython/AmberTools/Gaussian環境を有効にした後：

```bash
cd PMMA_TPS_OTf_MD
python scripts/server.py probe
cp config/server.example.json server.local.json
```

ここで`server.local.json`のCPU・メモリ・queue等を実環境に合わせて記入します。その後：

```bash
python scripts/server.py render --profile server.local.json --esp-method mk
qsub runs/qm/OTF_c01/job.pbs
```

**最初はOTF_c01だけ**を投入します。実行条件、正常終了、ESP点数・原子順・単位が確認できてから、残りを進めます。全ジョブを同時投入するスクリプトは用意していません。

## 次のファイルがまだ存在しない理由

採用済みの`PMMA19.itp`、`TPS.itp`、`OTF.itp`、`forcefield.itp`、`atomtypes.itp`は、実RESPとGAFF2変換・検査後に初めて生成します。完成した`.top`と混合系`init.gro`もその後です。今ある単分子`.gro`と`topol.top.in`だけではMDを始められません。

これはA案の完成を目指すための準備であり、A案の実計算・検証が完了した状態ではありません。
