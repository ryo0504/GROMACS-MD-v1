"""Pack only adopted, reviewed real force-field components; no MD execution."""
import argparse
import shutil
import math
from common import ROOT, read_json, write_json, execute, sha256


def pack(gmx):
    source=ROOT/'runs/adopted_input'
    adoption=read_json(source/'adoption.json')
    for name,digest in adoption['file_hashes'].items():
        if sha256(source/name)!=digest:raise ValueError('Adopted component changed: '+name)
    folder=ROOT/'runs/initial_box';folder.mkdir(parents=True,exist_ok=False)
    for path in source.iterdir():
        if path.suffix in ['.gro','.itp','.top']:shutil.copy2(path,folder/path.name)
    config=read_json(ROOT/'config/system.json');commands=[]
    expected=0;previous=None
    for name,output,seed in zip(['PMMA19','TPS','OTF'],['init0.gro','init1.gro','init.gro'],config['insert_seeds']):
        args=[gmx,'insert-molecules','-ci',name+'.gro','-nmol',str(config['counts'][name]),
              '-try','10000','-seed',str(seed),'-o',output]
        args+=['-f',previous] if previous else ['-box',*[str(x) for x in config['initial_box_nm']]]
        commands.append(execute(args,folder,folder/(name+'_insert.log')))
        n=int((folder/(name+'.gro')).read_text().splitlines()[1])
        expected+=n*config['counts'][name]
        lines=(folder/output).read_text().splitlines()
        if int(lines[1])!=expected or len(lines)!=expected+3:raise ValueError('Insertion count mismatch')
        previous=output
    if expected!=config['expected_atoms']:raise ValueError('Mixed system atom count differs from v2')
    write_json(folder/'packing_report.json',dict(status='PACKED_NOT_MINIMIZED',atoms=expected,
        total_charge_e=adoption['total_charge_e'],
        initial_density_g_cm3=adoption['total_mass_amu']*1.66053906660e-3/math.prod(config['initial_box_nm']),
        commands=commands,init_sha256=sha256(folder/'init.gro'),
        next='grompp/0-step preflight and opt; use the lab scheduler, not the login node'))
    print('Initial mixed box prepared:',folder/'init.gro')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--gmx',default='gmx')
    pack(p.parse_args().gmx)
