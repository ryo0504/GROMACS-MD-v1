# 未準備のRESP入力

ここは外部QM/ESP/RESP入力と電荷の出典を置くための未登録directoryです。テンプレートにはRESP電荷を同梱していません。2026-10-21に追加した[Notebook第5節のPySCFルート](../../../../examples/ethanol_nvt_comparison/README.md#このmacで第5節から進むpyscfルート)の実計算結果は`runs/ethanol27_pyscf/`に保存し、run-local configとhashで追跡します。

必要な記録：QM backend/version、method、basis、geometry、charge=0、multiplicity=1、ESP生成条件、RESP restraint/equivalence/stage条件、9原子の順序と電荷、元ファイルのhash。

単分子ASE原子順は `runs/ethanol27_initial/build_report.json` で確認します。電荷はその順序の9実数で `source_kind: charges` として読み込みます。AM1-BCCや仮の電荷に置き換えません。
