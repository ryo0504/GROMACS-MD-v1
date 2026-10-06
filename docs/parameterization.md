# GAFF / RESP

対象は有機低分子、簡単なoligomer/polymerモデルと有機混合物です。分子ごとのparameterを検証してからsystemを組み立てます。
力場familyは `gaff` または `gaff2` を明示します。charge modelはv1では `resp` のみです。

## 外部結果の受け入れ

`resp.source_file` の相対パスはsystem.yamlの位置から解決します。

| source_kind | 入力 | 続きの処理 |
|---|---|---|
| charges | atom順の空白区切りfitted charges | count/total charge確認後、antechamber -c rc |
| esp | Amber RESP形式ESP | 利用者指定fit_inputsを順にrespで実行 |
| gaussian_output | espgenで読めるGaussian ESP output | espgen → fit_inputsによるRESP |

Gaussianのすべての出力が対応するわけではありません。他QM backendの出力は外部でRESP ESP形式へ変換して `esp` として渡してください。
chargesを渡す場合もQM元データ・fitting条件をconfigとGitに残してください。
QM method/basis/backend、net charge、multiplicity、conformation、ESP手法を利用者が決めます。
`resp.atom_order_confirmed: true` は構造/QM/ESP/chargeのatom順とexplicit Hを確認してから設定します。
chargesの一致をatom数とtotal chargeだけから証明することはできません。

```yaml
resp:
  qm_backend: external
  method: HF          # 実際に使った条件を記録
  basis: 6-31G*
  source_kind: esp
  source_file: parameters/qm.esp
  atom_order_confirmed: true
  fit_inputs:
    - parameters/resp-stage1.in
    - parameters/resp-stage2.in
```

各RESP inputにはcharge、atom数、equivalence/restraint、stage2ならiqoptなどを科学的に適切に設定してください。
frameworkはconstraintsを自動生成しません。二段目以降に前段qoutを `-q` で渡します。
source不足時は明確なerrorとなり、chargeの別手法への変更は行いません。

## 検証

antechamberは保存済chargeを読み込み、parmchk2は選んだfamilyを使います。
frcmodに `ATTN` / `need revision` があれば停止します。自動推定された値は警告がなくても科学的レビューが必要です。
tleapはatom typing/parameterを組み立て、ParmEdでatom数、total charge、residue、bond/angle typeを確認します。
ACPYPEは既存AMBER topologyの変換モードを使い、chargeを再計算しません。
GROMACS前にはexpanded molecule count、atom/residue order、include、boxを確認し、最終的なmissing parameter/constraint検証はgromppで行います。
自動検証は力場の適用可能性やMDの科学的妥当性を証明しません。
