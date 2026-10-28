"""Single-point, matched-coordinate Amber/GROMACS conversion checks (server only)."""
import argparse
import re
import subprocess
from pathlib import Path
import numpy as np
from common import ROOT, read_json, write_json, execute, sha256


AMBER_INPUT='''Isolated component energy; not MD
&cntrl
 imin=1, maxcyc=0, ntb=0, igb=0, cut=999.0,
 ntc=1, ntf=1, ntpr=1,
/
'''
GROMACS_INPUT='''; Isolated-in-a-large-box single-point conversion check; not production.
integrator = md
nsteps = 0
dt = 0.001
pbc = xyz
cutoff-scheme = Verlet
nstlist = 1
verlet-buffer-tolerance = -1
rlist = 10.0
coulombtype = Cut-off
coulomb-modifier = None
rcoulomb = 10.0
vdwtype = Cut-off
vdw-modifier = None
rvdw = 10.0
DispCorr = no
constraints = none
tcoupl = no
pcoupl = no
gen-vel = no
continuation = no
comm-mode = Linear
nstcalcenergy = 1
nstenergy = 1
nstlog = 1
'''


def compare(name,gmx='gmx',case=None):
    import pyedr
    folder=case.resolve() if case else ROOT/'runs/forcefield'/name
    gro_name='coordinates.gro' if case else name+'.gro'
    if (folder/'energy_comparison.json').exists():raise ValueError('Energy comparison exists; inspect instead of overwriting')
    (folder/'amber_energy.in').write_text(AMBER_INPUT)
    (folder/'gmx_energy.mdp').write_text(GROMACS_INPUT)
    execute(['sander','-O','-i','amber_energy.in','-p','system.prmtop','-c','same_coordinates.inpcrd',
             '-o','amber_energy.out','-r','amber_energy.rst7'],folder,folder/'sander.command.log')
    execute([gmx,'grompp','-f','gmx_energy.mdp','-c',gro_name,'-p','isolated.top',
             '-o','energy.tpr','-po','energy_effective.mdp','-pp','energy_processed.top'],folder,folder/'energy_grompp.log')
    version=subprocess.run([gmx,'--version'],capture_output=True,text=True)
    args=[gmx,'mdrun','-deffnm','energy','-ntomp','1']
    if 'thread_mpi' in (version.stdout+version.stderr).lower():args+=['-ntmpi','1']
    execute(args,folder,folder/'energy_mdrun.log')
    edr=pyedr.edr_to_dict(str(folder/'energy.edr'))
    text=(folder/'amber_energy.out').read_text()
    def value(term):
        found=re.findall(re.escape(term)+r'\s*=\s*([-+0-9.EeDd]+)',text)
        if not found:raise ValueError('Amber energy missing: '+term)
        return float(found[-1].replace('D','E'))*4.184
    def gv(key):
        if key not in edr:raise ValueError('GROMACS energy missing: '+key)
        return float(edr[key][-1])
    dihkeys=[k for k in edr if 'Dih.' in k or 'Ryckaert' in k]
    pairs=[('bond',value('BOND'),gv('Bond')),('angle',value('ANGLE'),gv('Angle')),
           ('dihedral',value('DIHED'),sum(float(edr[k][-1]) for k in dihkeys)),
           ('LJ14',value('1-4 NB'),gv('LJ-14')),('Coulomb14',value('1-4 EEL'),gv('Coulomb-14')),
           ('LJ',value('VDWAALS'),gv('LJ (SR)')),('Coulomb',value('EELEC'),gv('Coulomb (SR)'))]
    rows=[]
    for term,a,g in pairs:
        tolerance=max(.01,1e-4*abs(a));diff=g-a
        rows.append(dict(term=term,amber_kJ_mol=a,gromacs_kJ_mol=g,
                         difference_kJ_mol=diff,tolerance_kJ_mol=tolerance,passed=abs(diff)<=tolerance))
    report=dict(passed=all(r['passed'] for r in rows),components=rows,
        top_sha256=sha256(folder/'isolated.top'),gro_sha256=sha256(folder/gro_name),
        prmtop_sha256=sha256(folder/'system.prmtop'),
        note='This is a single-geometry conversion check, not validation against experiment. Add conformers and ion-pair distances before mixed-system MD.')
    write_json(folder/'energy_comparison.json',report)
    print('Component energy comparison:',report['passed'],folder/'energy_comparison.json')
    if not report['passed']:raise ValueError('Resolve individual component differences; do not rely on total-energy cancellation')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('molecule',nargs='?',choices=['PMMA19','TPS','OTF'])
    p.add_argument('--case',type=Path);p.add_argument('--gmx',default='gmx');a=p.parse_args()
    if not a.molecule and not a.case:p.error('Supply molecule or --case')
    compare(a.molecule,a.gmx,a.case)
