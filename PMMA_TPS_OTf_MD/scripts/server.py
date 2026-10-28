"""Server discovery and explicit rendering. Does not submit any jobs."""
import argparse
import importlib.util
import importlib
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import warnings
from common import ROOT, read_json, write_json, sha256


def probe(output):
    tools={x:shutil.which(x) for x in ['g16','formchk','cubegen','espgen','resp','antechamber',
                                      'parmchk2','tleap','sander','acpype','gmx','gmx_mpi','qsub','qstat']}
    packages={}
    for name in ['numpy','scipy','rdkit','parmed','pyedr','matplotlib','pytest']:
        found=importlib.util.find_spec(name) is not None
        record={'found':found,'import_ok':False}
        if found:
            try:
                with warnings.catch_warnings(record=True) as seen:
                    warnings.simplefilter('always')
                    module=importlib.import_module(name)
                record.update(import_ok=True,version=getattr(module,'__version__',None),
                              warnings=[str(w.message) for w in seen])
            except Exception as e:record['error']=type(e).__name__+': '+str(e)
        packages[name]=record
    report=dict(platform=platform.platform(),machine=platform.machine(),cpu_count=os.cpu_count(),
                python=os.sys.executable,tools=tools,python_modules=packages,
                note='Read-only discovery; no Gaussian/QM/RESP/MD run and no job submission.')
    for executable in ['gmx','gmx_mpi']:
        if tools[executable]:
            cp=subprocess.run([tools[executable],'--version'],capture_output=True,text=True)
            report[executable+'_version']=cp.stdout+cp.stderr
    write_json(output,report)
    print(output)


def render(profile_path, esp_method):
    p=read_json(profile_path)
    if p.get('configured') is not True:
        raise ValueError('Set configured=true only after checking the lab server resource policy.')
    for key in ['nproc','memory_gb','pbs_memory_gb']:
        if type(p.get(key)) is not int or p[key]<=0: raise ValueError(f'Positive integer required: {key}')
    if p['pbs_memory_gb'] <= p['memory_gb']:
        raise ValueError('PBS memory must exceed Gaussian %mem (process/OS overhead).')
    if not re.fullmatch(r'[A-Za-z0-9_.@-]+',p.get('queue') or ''):
        raise ValueError('Set a valid lab PBS queue.')
    if not re.fullmatch(r'\d{1,3}:[0-5]\d:[0-5]\d',p.get('walltime') or ''):
        raise ValueError('Set PBS walltime as HH:MM:SS.')
    scratch=Path(p.get('scratch_root') or '')
    if not scratch.is_absolute() or not scratch.is_dir():
        raise ValueError('scratch_root must be an existing absolute directory on the server.')
    method=esp_method or p.get('esp_method')
    if method != 'mk':
        raise ValueError('The approved preparation protocol uses Gaussian16 MK ESP. Other grids require a separate protocol.')
    p['esp_method']=method
    target=ROOT/'runs/qm'
    if target.exists(): raise ValueError('runs/qm already exists; preserve results and choose a new run directory.')
    target.mkdir(parents=True)
    manifest=read_json(ROOT/'structures/manifest.json')
    q=shlex.quote
    python_cmd=shlex.join(shlex.split(p.get('python_command','python3')))
    g16_cmd=shlex.join(shlex.split(p.get('g16_command','g16')))
    root=str(ROOT.resolve())
    for job in manifest['jobs']:
        folder=target/job['id'];folder.mkdir()
        inp=(ROOT/'gaussian'/job['id']/'opt.gjf.in').read_text()
        inp=inp.replace('@NPROC@',str(p['nproc'])).replace('@MEM_GB@',str(p['memory_gb']))
        (folder/'opt.gjf').write_text(inp)
        if method=='mk':
            # Log ESP is deliberately requested for AmberTools espgen.
            # The grid differs from the v2 Fibonacci grid; record it as a protocol variant.
            esp=(ROOT/'gaussian'/job['id']/'esp.gjf.in').read_text()
            (folder/'esp.gjf').write_text(esp.replace('@NPROC@',str(p['nproc'])).replace('@MEM_GB@',str(p['memory_gb'])))
        setup=''
        for key in ['gaussian_setup','amber_setup']:
            if p.get(key):
                setup+='source '+q(str(p[key]))+'\n'
        script=(f"#!/usr/bin/env bash\n#PBS -N {job['id']}\n#PBS -q {p['queue']}\n"
                f"#PBS -l select=1:ncpus={p['nproc']}:ompthreads={p['nproc']}:mem={p['pbs_memory_gb']}gb\n"
                f"#PBS -l walltime={p['walltime']}\n#PBS -j oe\nset -euo pipefail\n"
                +setup+f'cd {q(str(folder.resolve()))}\n'
                +f'export GAUSS_SCRDIR="$(mktemp -d {q(str(scratch / (job["id"]+".XXXXXX")))})"\n'
                +'trap \'status=$?; if [ "$status" -eq 0 ]; then rm -rf -- "$GAUSS_SCRDIR"; else printf "Failed-job scratch retained: %s\\n" "$GAUSS_SCRDIR" >&2; fi\' EXIT\n'
                +'test ! -e opt.log || { echo "opt.log exists; inspect before restarting" >&2; exit 1; }\n'
                +f'{g16_cmd} < opt.gjf > opt.log\n'
                +f'{python_cmd} {q(root+"/scripts/gaussian_results.py")} check-log opt.log --optimization\n')
        if method=='mk':
            script+=f'{g16_cmd} < esp.gjf > esp.log\n'
            script+=f'{python_cmd} {q(root+"/scripts/gaussian_results.py")} check-log esp.log\n'
        script+=f'formchk optimized.chk optimized.fchk\n'
        script+=f'{python_cmd} {q(root+"/scripts/gaussian_results.py")} extract {q(job["id"])} --method {method}\n'
        (folder/'job.pbs').write_text(script)
        write_json(folder/'job.json',dict(**job,esp_method=method,profile_sha256=sha256(profile_path),
            geometry_input_sha256=sha256(folder/'opt.gjf'),status='PREPARED_NOT_SUBMITTED'))
    write_json(ROOT/'runs/rendered_profile.json',p)
    write_json(ROOT/'reports/gaussian_protocol_selection.json',dict(esp_method=method,
        v2_difference='Gaussian16 MK grid instead of the v2 Fibonacci grid (user-approved 2026-10-28)',
        status='PREPARED_NOT_RUN'))
    print('Rendered jobs:',len(manifest['jobs']))
    print('First job (submit on compute server only after review):',target/'OTF_c01/job.pbs')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='action',required=True)
    a=sub.add_parser('probe');a.add_argument('--output',type=Path,default=ROOT/'runs/server_probe.json')
    a=sub.add_parser('render');a.add_argument('--profile',type=Path,required=True)
    a.add_argument('--esp-method',choices=['mk'])
    args=parser.parse_args()
    if args.action=='probe':probe(args.output)
    else:render(args.profile,args.esp_method)


if __name__=='__main__':main()
