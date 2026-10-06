"""Register this host's Workbench kernel, including external tool discovery."""
from pathlib import Path
import json
import os
import subprocess
import sys
from jupyter_client.kernelspec import KernelSpecManager

subprocess.run([sys.executable, '-m', 'ipykernel', 'install', '--user', '--name', 'gmx-md',
                '--display-name', 'Python (GROMACS MD Workbench)'], check=True)
path = Path(KernelSpecManager().get_kernel_spec('gmx-md').resource_dir) / 'kernel.json'
data = json.loads(path.read_text())
data['env'] = {
    'PATH': os.environ['PATH'],
    'CONDA_PREFIX': sys.prefix,
    'CONDA_DEFAULT_ENV': 'gmx-md',
    'AMBERHOME': os.environ.get('AMBERHOME', sys.prefix),
}
for key in ('MAMBA_ROOT_PREFIX', 'GMX_BIN'):
    if os.environ.get(key):
        data['env'][key] = os.environ[key]
path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
print('Kernel:', path)
print('Python:', data['argv'][0])
print('External GROMACS:', data['env'].get('GMX_BIN', 'PATH lookup'))
