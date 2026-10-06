# Ethanol 27 benchmark scientific inputs

教材と実行手順は [examples/ethanol_nvt_comparison/README.md](../../examples/ethanol_nvt_comparison/README.md)。

`system.yaml` がこの教材の科学条件です。RESPや非結合条件のnullは未決定です。未入力ではparameterize/prepareが停止します。buildとRDF演習は独立して実行できます。

`parameters/resp/` は外部QM/ESP/RESP入力・電荷・原子順序の確認記録を置く場所です。大きいQM outputはruns等に置き、configとhashで追跡してください。`analysis/reference/` は将来のMACE比較元の出典・step・原子mapping等を記録する場所です。現在比較データはありません。

実行出力はrepository rootの `runs/ethanol27_*` に保存し、Git管理対象にしません。
