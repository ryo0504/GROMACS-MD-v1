# Cluster

`detected/` にprobe結果、`slurm/` に資源未指定のテンプレートを配置します。
reportを読み、既存module、GPU/MPI build、scheduler、filesystemとsite policyから設定してください。
jobはstageごとに実行します。production templateはNPT完了済みrunを前提とします。
テンプレートは未設定状態でsubmitできません。ログ出力directoryはsbatch前に作成してください。
詳細: [server setup](../docs/server_setup.md)。
