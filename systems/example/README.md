# Example: ethanol

CCOをexplicit H付き3D構造にします。system.yamlの条件は教育用の例で、RESP sourceとbulk systemは未提供です。
`parameterization_file` はQMと同じatom順/geometryの構造に置き換えます。
`parameters/` にレビュー済みRESP charges/input、frcmod、MOL2、top/itpなどを保存できます。
SMILES→3D生成のみでbulk ethanolの科学的モデルは完成しません。
packing/分子数/box/density、適用力場、QM/RESP、MDPは利用者の科学的判断で確定してください。
