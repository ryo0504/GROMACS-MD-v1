"""Generate explicit structures and input files; no QM, RESP or MD execution."""
import argparse
import csv
from collections import defaultdict
import importlib.metadata
from pathlib import Path
import sys

from common import ROOT, write_json, sha256


ROLES = ['backbone_head_C', 'backbone_substituted_C', 'alpha_methyl_C',
         'carbonyl_C', 'carbonyl_O', 'ester_O', 'methoxy_C']


def polymer(n):
    from rdkit import Chem
    import random
    rw = Chem.RWMol()
    previous = None
    for unit in range(n):
        monomer = Chem.MolFromSmiles('CC(C)(C(=O)OC)')
        offset = rw.GetNumAtoms()
        rw.InsertMol(monomer)
        for j, role in enumerate(ROLES):
            atom = rw.GetAtomWithIdx(offset+j)
            atom.SetIntProp('unit', unit+1)
            atom.SetProp('role', role)
        if previous is not None:
            rw.AddBond(previous, offset, Chem.BondType.SINGLE)
        previous = offset+1
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    centers = [i for i, _ in Chem.FindMolChiralCenters(mol, includeUnassigned=True)]
    rng = random.Random(42)
    for i in centers:
        mol.GetAtomWithIdx(i).SetChiralTag(rng.choice([
            Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW]))
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            parent = atom.GetNeighbors()[0]
            hs = sorted(a.GetIdx() for a in parent.GetNeighbors() if a.GetSymbol() == 'H')
            atom.SetIntProp('unit', parent.GetIntProp('unit'))
            atom.SetProp('role', parent.GetProp('role') + '_H' + str(hs.index(atom.GetIdx())+1))
    return mol


def mapping(mol, name, degree=None):
    from rdkit import Chem
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True))
    rows = []
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        # 5 columns in GRO atom-name field, stable names across formats.
        atom_name = atom.GetSymbol().upper() + str(i+1)
        methyl = atom.GetSymbol() == 'C' and sum(x.GetSymbol() == 'H' for x in atom.GetNeighbors()) == 3
        parent = atom.GetNeighbors()[0] if atom.GetSymbol() == 'H' else atom
        parent_methyl = parent.GetSymbol() == 'C' and sum(x.GetSymbol() == 'H' for x in parent.GetNeighbors()) == 3
        if degree:
            unit = atom.GetIntProp('unit')
            group = 'head' if unit == 1 else ('tail' if unit == degree else 'internal')
            role = atom.GetProp('role')
            eq_role = role.rsplit('_H', 1)[0] + '_H' if atom.GetSymbol() == 'H' and parent_methyl else role
            key = group + ':' + eq_role
        else:
            unit, group = 1, name
            role = atom.GetSymbol() + str(i+1)
            key = name + ':sym' + str(ranks[i])
            if name == 'OTF' and atom.GetSymbol() in ('O', 'F'):
                key = 'OTF:equivalent_' + atom.GetSymbol()
        atom.SetAtomMapNum(i+1)
        rows.append(dict(index=i+1, element=atom.GetSymbol(), atomic_number=atom.GetAtomicNum(),
                         atom_name=atom_name, unit=unit, template=group, role=role,
                         equivalence_key=key, stage2_free=bool(methyl or (atom.GetSymbol() == 'H' and parent_methyl)),
                         formal_charge=atom.GetFormalCharge(),
                         cip=atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None,
                         neighbors=[a.GetIdx()+1 for a in atom.GetNeighbors()]))
    return rows


def write_structure(mol, conf, path, rows):
    from rdkit import Chem
    mol = Chem.Mol(mol, confId=conf)
    path.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(path.with_suffix('.sdf'))) as writer:
        writer.write(mol)
    xyz = mol.GetConformer().GetPositions()
    name = path.stem.split('_')[0]
    resname = 'PMM' if name.startswith('PMMA') else name
    path.with_suffix('.xyz').write_text(str(len(rows)) + '\nInitial RDKit conformer; NOT QM optimized\n' +
        '\n'.join(f"{r['element']} {x:.10f} {y:.10f} {z:.10f}" for r, (x,y,z) in zip(rows, xyz)) + '\n')
    gro = ['Initial structure only; force field and charges NOT_GENERATED', str(len(rows))]
    for r, (x,y,z) in zip(rows, xyz):
        gro.append(f"{1:5d}{resname:<5}{r['atom_name']:>5}{r['index']:5d}{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}")
    path.with_suffix('.gro').write_text('\n'.join(gro) + '\n  10.00000  10.00000  10.00000\n')


def make_gaussian_templates(job, rows, xyz):
    """Generate route-neutral optimization; ESP route is selected on the server."""
    folder = ROOT/'gaussian'/job['id']
    folder.mkdir(parents=True, exist_ok=True)
    header = '%nprocshared=@NPROC@\n%mem=@MEM_GB@GB\n%chk=optimized.chk\n'
    route = '#p HF/6-31G* Opt=(Tight,MaxCycles=100) SCF=(Conver=10,MaxCycle=100) NoSymm\n\n'
    title = job['id'] + ': RHF/6-31G* geometry optimization; atom order fixed\n\n'
    geometry = f"{job['charge']} 1\n" + '\n'.join(
        f"{r['element']} {x:.10f} {y:.10f} {z:.10f}" for r,(x,y,z) in zip(rows,xyz)) + '\n\n'
    (folder/'opt.gjf.in').write_text(header + route + title + geometry)
    (folder/'esp.gjf.in').write_text(header +
        '#p HF/6-31G* Geom=AllCheck Guess=Read SCF=(Conver=10,MaxCycle=100)\n'
        ' NoSymm Pop=MK IOp(6/33=2,6/41=4,6/42=6)\n\n')


def make_mdp():
    common = {'integrator':'md','dt':'0.001','tinit':0,'pbc':'xyz','cutoff-scheme':'Verlet',
        'nstlist':20,'verlet-buffer-tolerance':.005,'coulombtype':'PME','rcoulomb':1.0,
        'ewald-rtol':'1e-5','pme-order':4,'fourierspacing':.12,'vdwtype':'Cut-off','rvdw':1.0,
        'vdw-modifier':'Potential-shift','DispCorr':'EnerPres','constraints':'h-bonds',
        'constraint-algorithm':'lincs','lincs-order':4,'lincs-iter':1,'tcoupl':'v-rescale',
        'tc-grps':'System','tau-t':1.0,'pcoupl':'c-rescale','pcoupltype':'isotropic',
        'tau-p':5.0,'compressibility':'4.5e-5','comm-mode':'Linear','comm-grps':'System',
        'nstcomm':100,'nstcalcenergy':100,'nstenergy':1000,'nstlog':1000,
        'nstxout':0,'nstvout':0,'nstfout':0,'nstxout-compressed':1000,
        'compressed-x-precision':10000,'gen-vel':'no','continuation':'yes','ld-seed':201,'annealing':'no'}
    stages = [('heat',200000,300,None,'opt'),('eq_1',500000,1400,100,'heat'),
              ('eq_2',500000,1400,10,'eq_1'),('eq_3',5000000,1400,1,'eq_2'),
              ('prod',20000000,1400,1,'eq_3'),('cool',11000000,1400,1,'prod'),
              ('eq_rt',20000000,300,1,'cool'),('prod_rt',20000000,300,1,'eq_rt'),
              ('nvt_rt',200000,300,None,'eq_rt'),('prod_rt_nvt',20000000,300,None,'nvt_rt')]
    ends = {'opt':0}; schedule = {}
    folder = ROOT/'mdp'; folder.mkdir(exist_ok=True)
    em = {k:common[k] for k in ['pbc','cutoff-scheme','nstlist','coulombtype','rcoulomb',
                               'vdwtype','rvdw','vdw-modifier','DispCorr','ewald-rtol','pme-order','fourierspacing']}
    em.update(integrator='steep',nsteps=50000,emtol=100,emstep=.01,constraints='none')
    def output(name, values):
        (folder/(name+'.mdp')).write_text('; Prepared input. QM/RESP/topology validation still required.\n' +
                                          '\n'.join(f'{k:<26} = {v}' for k,v in values.items())+'\n')
    output('opt',em)
    pre=dict(common)
    for key in ['tau-t','tc-grps','pcoupltype','tau-p','compressibility']:
        pre.pop(key,None)
    pre.update({'nsteps':0,'init-step':0,'tcoupl':'no','pcoupl':'no',
                'gen-vel':'no','continuation':'no','constraints':'none',
                'nstxout':1,'nstfout':1,'nstxout-compressed':0,'nstcalcenergy':1,'nstenergy':1,'nstlog':1})
    output('preflight',pre)
    for name,nsteps,temp,pressure,parent in stages:
        d = dict(common); start = ends[parent]
        d.update({'init-step':start,'nsteps':nsteps,'ref-t':temp})
        if pressure is None:
            d['pcoupl']='no'
            for k in ['pcoupltype','tau-p','compressibility']: d.pop(k)
        else: d['ref-p']=pressure
        if name == 'heat':
            d.update({'gen-vel':'yes','continuation':'no','gen-temp':300,'gen-seed':1,
                      'annealing':'single','annealing-npoints':2,'annealing-time':'0 200',
                      'annealing-temp':'300 1400'})
        if name == 'cool':
            times = [0,2000,4000,6000,8000,9000,9500,10000,10500,11000]
            d.update({'annealing':'single','annealing-npoints':10,
                      'annealing-time':' '.join(f'{start*.001+t:g}' for t in times),
                      'annealing-temp':'1400 1200 1000 800 600 500 450 400 350 300'})
        output(name,d)
        ends[name]=start+nsteps
        schedule[name]=dict(parent=parent,init_step=start,nsteps=nsteps,dt_ps=.001,
                            start_ps=start*.001,end_ps=ends[name]*.001,temperature_K=temp,
                            pressure_bar=pressure,optional=name in ['nvt_rt','prod_rt_nvt'])
    write_json(ROOT/'config/protocol.json',dict(stages=schedule,
        main_duration_ns=77.2,with_rt_nvt_branch_duration_ns=97.4,
        note='If any upstream stage is extended, regenerate init-step/annealing-time; do not reuse stale MDP.'))


def prepare():
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    if (ROOT/'structures/manifest.json').exists():
        raise SystemExit('Already generated. Preserve inputs; use a new directory for regeneration.')
    for f in ['config','structures','gaussian','topology','reports','templates/pbs','tests']:
        (ROOT/f).mkdir(parents=True, exist_ok=True)
    jobs=[]; models={}; sequences={}
    for name,degree,count,charge in [('PMMA19',19,1,0),('PMMA3',3,3,0),('PMMA5',5,3,0),
                                      ('TPS',None,3,1),('OTF',None,1,-1)]:
        mol=polymer(degree) if degree else Chem.AddHs(Chem.MolFromSmiles(
            'c1ccc([S+](c2ccccc2)c2ccccc2)cc1' if name=='TPS' else 'O=S(=O)([O-])C(F)(F)F'))
        rows=mapping(mol,name,degree)
        params=AllChem.ETKDGv3();params.randomSeed=43;params.numThreads=1;params.enforceChirality=True
        params.pruneRmsThresh=-1
        if degree == 19:
            # Bounds-matrix embedding of a highly stereospecified chain often fails.
            # Random-coordinate distance geometry retains the same explicit stereo constraints.
            params.useRandomCoords=True
        params.timeout=30
        ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=count,params=params))
        if len(ids)!=count: raise ValueError(f'Embedding failed: {name}')
        optim=[]
        for conf in ids:
            if degree and AllChem.MMFFHasAllMoleculeParams(mol):
                status=AllChem.MMFFOptimizeMolecule(mol,confId=conf,maxIters=300)
                optim.append({'method':'MMFF94 (initial coordinates only)','status':status})
            else: optim.append({'method':'ETKDG only; not energy optimized','status':None})
        for r in rows:
            assert len(r['atom_name']) <= 5
        write_json(ROOT/'structures'/(name+'_atom_map.json'),rows)
        with (ROOT/'structures'/(name+'_atom_map.csv')).open('w',newline='') as f:
            fields=[k for k in rows[0] if k!='neighbors']
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
            w.writerows({k:r[k] for k in fields} for r in rows)
        for j,conf in enumerate(ids,1):
            job_id=f'{name}_c{j:02d}'
            path=ROOT/'structures'/job_id
            write_structure(mol,conf,path,rows)
            if name!='PMMA19':
                purpose='evaluation' if name in ['PMMA3','PMMA5'] and j==3 else 'fit'
                job=dict(id=job_id,molecule=name,charge=charge,multiplicity=1,
                         atom_count=len(rows),purpose=purpose,
                         sdf=str(path.with_suffix('.sdf').relative_to(ROOT)))
                make_gaussian_templates(job,rows,mol.GetConformer(conf).GetPositions())
                jobs.append(job)
        models[name]=dict(degree=degree,formula=rdMolDescriptors.CalcMolFormula(mol),
                          atoms=mol.GetNumAtoms(),formal_charge=Chem.GetFormalCharge(mol),
                          initial_conformers=optim,
                          stereo=[{'index':r['index'],'unit':r['unit'],'cip':r['cip']} for r in rows if r['cip']])
    write_json(ROOT/'structures/manifest.json',dict(models=models,jobs=jobs,
        rdkit_version=importlib.metadata.version('rdkit'),
        status='INITIAL_STRUCTURES_ONLY; NO_QM_OR_RESP',seed=43,
        note='PMMA3/5/19 each have a fixed seeded sequence; they are not claimed to be identical tacticity environments.'))
    write_json(ROOT/'config/system.json',dict(schema_version=2,polymer='PMMA19',
        counts={'PMMA19':100,'TPS':83,'OTF':83},initial_box_nm=[10.,10.,10.],
        force_field='gaff2',charge_model='RESP',charge_scale=1.,
        net_charge=0,expected_atoms=32186,insert_seeds=[1,1,1],
        forcefield_status='NOT_GENERATED',qm_backend='Gaussian16',esp_method='mk'))
    write_json(ROOT/'config/server.example.json',dict(configured=False,nproc=None,memory_gb=None,
        pbs_memory_gb=None,queue=None,walltime=None,gaussian_setup=None,amber_setup=None,
        python_command='python3',g16_command='g16',gmx_command='gmx',
        scratch_root=None,esp_method='mk'))
    (ROOT/'topology/topol.top.in').write_text('; TEMPLATE: includes are generated only after real RESP/GAFF2 validation.\n'
        '#include "forcefield.itp"\n#include "atomtypes.itp"\n#include "PMMA19.itp"\n'
        '#include "TPS.itp"\n#include "OTF.itp"\n\n[ system ]\nPMMA19 TPS triflate v2\n\n'
        '[ molecules ]\nPMMA19 100\nTPS 83\nOTF 83\n')
    make_mdp()
    spec=ROOT/'references/PMMA_TPS_OTf_MD_tutorial_spec_v2.md'
    if not spec.exists():spec=ROOT.parent/'MD習熟/PMMA_TPS_OTf_MD_tutorial_spec_v2.md'
    write_json(ROOT/'reports/preparation_status.json',dict(date='2026-10-28',
        prepared=['canonical_initial_structures','atom_maps','Gaussian_optimization_templates','MDP','configuration'],
        not_run=['Gaussian16','ESP','RESP','GAFF2','Amber_to_Gromacs_energy_check','packing','MD','PBS'],
        spec_sha256=sha256(spec)))
    print('Prepared initial structures, Gaussian templates and MDP. No QM/RESP/MD executed.')


if __name__=='__main__':
    prepare()
