# Architecture

科学的条件は YAML、MDP、RESP input、topology とコードで管理し、実行 binary はホストごとに選択します。
Python 3.11 を共通基盤とし、GROMACS は外部実行ファイルです。

| 層 | 責任 |
|---|---|
| notebooks | 探索・レビュー・図 |
| src/mdworkbench | 構造、parameterization、検証、解析、provenance |
| scripts/run, cluster/slurm | 実行・資源割当 |
| config, systems, mdp | 明示した科学的条件 |
| env | portable environment と platform 別 lock |
| runs (ignored) | 独立した run、入力 snapshot、結果 |

重い Python dependency は遅延 import。doctor と外部 command detection は分子ライブラリなしでも実行可能です。
RESP は外部 QM/ESP と利用者指定の fitting input、または atom 順序を確認済みの fitted charge を受け取ります。
QM 自動実行、溶媒/混合物自動 packing、複雑な polymer 構築は v1 の範囲外です。
run は topology の include を snapshot し、grompp の処理済み topology/MDP を残します。
grompp warning の無視や charge model の fallback は行いません。
