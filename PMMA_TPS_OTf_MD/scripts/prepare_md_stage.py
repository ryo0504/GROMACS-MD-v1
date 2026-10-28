"""Prepare one reviewed MD stage and its commands; never submit or run mdrun."""
import argparse
from pathlib import Path
import shlex
import shutil
import subprocess
import re
from common import ROOT, read_json, write_json, sha256


def prepare(stage,parent_gro,parent_cpt,gmx):
    conversion=read_json(ROOT/'reports/conversion_review.json')
    required=['single_molecules','multiple_polymer_geometries','ion_pair_distances_and_orientations']
    if conversion.get('approved') is not True or not all(conversion.get(k) is True for k in required):
        raise ValueError('Complete the v2 component-energy/parameter review matrix before preparing MD.')
    evidence=conversion.get('evidence_files',[])
    if len(evidence)<11 or len({e['path'] for e in evidence})!=len(evidence):
        raise ValueError('The 3 single-molecule and 8 additional energy cases need distinct real evidence.')
    for item in evidence:
        if sha256(ROOT/item['path'])!=item['sha256']:
            raise ValueError('Conversion evidence changed after approval')
    protocol=read_json(ROOT/'config/protocol.json')
    if stage!='opt' and stage not in protocol['stages']:raise ValueError('Unknown MD stage')
    source=ROOT/'runs/adopted_input'
    adoption=read_json(source/'adoption.json')
    for name,digest in adoption['file_hashes'].items():
        if sha256(source/name)!=digest:raise ValueError('Adopted input hash changed')
    if stage=='opt':
        parent_gro=ROOT/'runs/initial_box/init.gro'
        packed=read_json(ROOT/'runs/initial_box/packing_report.json')
        if sha256(parent_gro)!=packed['init_sha256']:raise ValueError('Packed initial box changed')
    if parent_gro is None or not parent_gro.is_file():raise ValueError('Supply the actual adopted parent GRO file')
    if stage not in ['opt','heat'] and (parent_cpt is None or not parent_cpt.is_file()):
        raise ValueError('Supply the actual parent checkpoint; do not regenerate velocities')
    if stage in ['opt','heat'] and parent_cpt is not None:
        raise ValueError('The first EM/MD stage does not import a prior checkpoint')
    dump=None
    if parent_cpt:
        cp=subprocess.run([gmx,'dump','-cp',str(parent_cpt.resolve())],capture_output=True,text=True)
        if cp.returncode:raise ValueError('Cannot read parent checkpoint with the selected GROMACS')
        dump=cp.stdout+cp.stderr
        found=re.search(r'(?m)^\s*step\s*=\s*(\d+)',dump)
        if not found or int(found.group(1))!=protocol['stages'][stage]['init_step']:
            raise ValueError('Parent checkpoint step differs from planned init-step. Regenerate the schedule/annealing-time after any extension.')
    folder=ROOT/'runs/md'/stage;folder.mkdir(parents=True,exist_ok=False)
    if dump:(folder/'parent_checkpoint_dump.txt').write_text(dump)
    for path in source.iterdir():
        if path.suffix in ['.itp','.top']:shutil.copy2(path,folder/path.name)
    shutil.copy2(ROOT/'mdp'/(stage+'.mdp'),folder/(stage+'.mdp'))
    q=shlex.quote
    prep=[gmx,'grompp','-f',stage+'.mdp','-c',str(parent_gro.resolve()),'-p','topol.top',
          '-o',stage+'.tpr','-po',stage+'_effective.mdp','-pp',stage+'_processed.top']
    if parent_cpt:prep+=['-t',str(parent_cpt.resolve())]
    run=[gmx,'mdrun','-deffnm',stage,'-cpt','5']
    if parent_cpt:run+=['-cpi',str(parent_cpt.resolve()),'-noappend']
    text=('#!/usr/bin/env bash\nset -euo pipefail\ncd '+q(str(folder.resolve()))+'\n'
          +'# Run only inside the reviewed PBS allocation and GROMACS/MPI environment.\n'
          +'# Match -maxh and MPI launcher to the successful lab job; this script does not select resources.\n'
          +shlex.join(prep)+'\n'+shlex.join(run)+'\n')
    (folder/'commands.sh').write_text(text)
    write_json(folder/'stage_preparation.json',dict(stage=stage,status='PREPARED_NOT_EXECUTED',
        parent_gro=str(parent_gro.resolve()),parent_cpt=str(parent_cpt.resolve()) if parent_cpt else None,
        parent_gro_sha256=sha256(parent_gro),mdp_sha256=sha256(folder/(stage+'.mdp')),
        expected=protocol['stages'].get(stage),
        reminder='Verify actual checkpoint step/time against init-step before running grompp; if upstream was extended, regenerate protocol.'))
    print('Prepared commands, not submitted:',folder/'commands.sh')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage')
    p.add_argument('--parent-gro',type=Path);p.add_argument('--parent-cpt',type=Path);p.add_argument('--gmx',default='gmx')
    a=p.parse_args();prepare(a.stage,a.parent_gro,a.parent_cpt,a.gmx)
