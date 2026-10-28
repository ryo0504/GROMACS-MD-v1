"""Check actual Gaussian results and extract real ESP; never synthesize results."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np
from common import ROOT, BOHR_A, read_json, write_json, read_esp, execute, sha256


def check_log(path, optimization=False):
    text=Path(path).read_text(errors='replace')
    if 'Normal termination of Gaussian' not in text[-4000:] or 'Error termination' in text:
        raise ValueError(f'Gaussian did not terminate normally: {path}')
    if optimization and not ('Optimization completed.' in text or 'Stationary point found.' in text):
        raise ValueError(f'Geometry optimization completion not found: {path}')
    if 'SCF Done:' not in text:raise ValueError('No converged SCF energy found')
    return text


def fchk_array(path,key):
    lines=Path(path).read_text().splitlines()
    for i,line in enumerate(lines):
        if line.startswith(key) and 'N=' in line:
            n=int(line.split('N=')[1]);values=[]
            for l in lines[i+1:]:
                values.extend(l.split())
                if len(values)>=n:return np.array([float(x.replace('D','E')) for x in values[:n]])
    raise ValueError(f'fchk array missing: {key}')


def fchk_integer(path,key):
    for line in Path(path).read_text().splitlines():
        if line.startswith(key):
            match=re.search(r'\bI\s+(-?\d+)\s*$',line)
            if match:return int(match.group(1))
    raise ValueError(f'fchk integer missing: {key}')


def write_esp(path,coords,grid,potential):
    if len(coords)>99999 or len(grid)>99999:raise ValueError('RESP 5-column dimensions exceeded')
    with Path(path).open('x') as h:
        h.write(f'{len(coords):5d}{len(grid):5d}{0:5d}\n')
        for x,y,z in coords:h.write(' '*17+f'{x:16.7E}{y:16.7E}{z:16.7E}\n')
        for v,(x,y,z) in zip(potential,grid):h.write(f' {v:16.7E}{x:16.7E}{y:16.7E}{z:16.7E}\n')


def extract(job_id,method):
    from rdkit import Chem
    manifest=read_json(ROOT/'structures/manifest.json')
    job=next(j for j in manifest['jobs'] if j['id']==job_id)
    folder=ROOT/'runs/qm'/job_id
    check_log(folder/'opt.log',True)
    meta=read_json(folder/'job.json')
    if meta['esp_method']!=method:raise ValueError('ESP method differs from rendered input')
    maps=read_json(ROOT/'structures'/(job['molecule']+'_atom_map.json'))
    numbers=fchk_array(folder/'optimized.fchk','Atomic numbers').astype(int)
    xyz_bohr=fchk_array(folder/'optimized.fchk','Current cartesian coordinates').reshape(-1,3)
    if fchk_integer(folder/'optimized.fchk','Charge')!=job['charge']:
        raise ValueError('Gaussian charge does not match the specified molecule')
    if fchk_integer(folder/'optimized.fchk','Multiplicity')!=1:
        raise ValueError('Gaussian multiplicity does not match the closed-shell singlet protocol')
    if numbers.tolist()!=[r['atomic_number'] for r in maps]:raise ValueError('Gaussian atomic number/order mismatch')
    mol=Chem.SDMolSupplier(str(ROOT/job['sdf']),removeHs=False)[0]
    expected_stereo={a.GetIdx():a.GetProp('_CIPCode') for a in mol.GetAtoms() if a.HasProp('_CIPCode')}
    xyz_A=xyz_bohr*BOHR_A
    for i,pos in enumerate(xyz_A):mol.GetConformer().SetAtomPosition(i,pos.tolist())
    # Use covalent radii to check intended bonds, without re-inventing charged sulfur bond orders.
    pt=Chem.GetPeriodicTable()
    for b in mol.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        radius=pt.GetRcovalent(int(numbers[i]))+pt.GetRcovalent(int(numbers[j]))
        distance=np.linalg.norm(xyz_A[i]-xyz_A[j])
        if not .6*radius<distance<1.5*radius:raise ValueError(f'Unexpected optimized bond length: {i+1}-{j+1}')
    Chem.AssignStereochemistryFrom3D(mol,replaceExistingTags=True)
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    for i,cip in expected_stereo.items():
        a=mol.GetAtomWithIdx(i)
        if not a.HasProp('_CIPCode') or a.GetProp('_CIPCode')!=cip:raise ValueError(f'Stereochemistry changed: atom {i+1}')
    if (folder/'result.json').exists():raise ValueError('Result already extracted; inspect before replacing.')
    esp=folder/'qm.esp'
    if method != 'mk':raise ValueError('This bundle is configured for the approved MK-grid protocol only')
    check_log(folder/'esp.log')
    execute(['espgen','-i','esp.log','-o','qm.esp'],folder,folder/'espgen.log')
    nuclei,grid,values=read_esp(esp)
    if nuclei.shape!=xyz_bohr.shape or not np.allclose(nuclei,xyz_bohr,atol=2e-4):
        raise ValueError('ESP/fchk geometry mismatch (orientation, units, order or wrong checkpoint)')
    with Chem.SDWriter(str(folder/'optimized.sdf')) as w:w.write(mol)
    write_json(folder/'result.json',dict(job_id=job_id,status='QM_AND_ESP_EXTRACTED_NOT_FORCEFIELD_VALIDATED',
        esp_method=method,atom_count=len(nuclei),esp_points=len(grid),esp_sha256=sha256(esp),
        input_sha256=sha256(folder/'opt.gjf'),fchk_sha256=sha256(folder/'optimized.fchk'),
        optimized_sdf_sha256=sha256(folder/'optimized.sdf')))
    print('Validated real QM/ESP:',job_id,len(grid),'points')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('check-log');a.add_argument('path',type=Path);a.add_argument('--optimization',action='store_true')
    a=s.add_parser('extract');a.add_argument('job');a.add_argument('--method',choices=['mk'],required=True)
    a=p.parse_args()
    if a.action=='check-log':check_log(a.path,a.optimization)
    else:extract(a.job,a.method)
