# 未準備のRESP入力

ここには実際のQM/ESP/RESP入力と電荷の出典を記録します。現時点でRESP電荷はありません。

必要な記録：QM backend/version、method、basis、geometry、charge=0、multiplicity=1、ESP生成条件、RESP restraint/equivalence/stage条件、9原子の順序と電荷、元ファイルのhash。

単分子ASE原子順は `runs/ethanol27_initial/build_report.json` で確認します。電荷はその順序の9実数で `source_kind: charges` として読み込みます。AM1-BCCや仮の電荷に置き換えません。
