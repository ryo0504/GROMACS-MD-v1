"""Lightweight preparation tests. Synthetic ESP fixtures are NOT molecular results."""
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import read_esp,read_json
from gaussian_results import write_esp,check_log
from fit_resp import controls,read_qout,maps_for


def test_structures_formula_stereo_and_mapping():
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    expected={'PMMA19':('C95H154O38',287,0),'PMMA3':('C15H26O6',47,0),
              'PMMA5':('C25H42O10',77,0),'TPS':('C18H15S+',34,1),'OTF':('CF3O3S-',8,-1)}
    manifest=read_json(ROOT/'structures/manifest.json')
    for name,(formula,n,q) in expected.items():
        mapping=read_json(ROOT/'structures'/(name+'_atom_map.json'))
        for path in sorted((ROOT/'structures').glob(name+'_c*.sdf')):
            mol=Chem.SDMolSupplier(str(path),removeHs=False)[0]
            assert mol is not None and mol.GetNumAtoms()==n
            assert rdMolDescriptors.CalcMolFormula(mol)==formula
            assert Chem.GetFormalCharge(mol)==q
            assert [a.GetAtomMapNum() for a in mol.GetAtoms()]==list(range(1,n+1))
            assert [a.GetSymbol() for a in mol.GetAtoms()]==[r['element'] for r in mapping]
            before={a.GetIdx():a.GetProp('_CIPCode') for a in mol.GetAtoms() if a.HasProp('_CIPCode')}
            Chem.AssignStereochemistryFrom3D(mol,replaceExistingTags=True)
            Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
            assert all(mol.GetAtomWithIdx(i).GetProp('_CIPCode')==cip for i,cip in before.items())
            xyz=mol.GetConformer().GetPositions();assert np.isfinite(xyz).all()
            for bond in mol.GetBonds():
                a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
                assert .6 < np.linalg.norm(xyz[a]-xyz[b]) < 2.2
            # Nonbonded initial contacts below 0.7 Å indicate an unusable starting conformer.
            d=np.linalg.norm(xyz[:,None]-xyz[None],axis=2);np.fill_diagonal(d,np.inf)
            assert d.min()>.7
        assert len(mapping)==n
        assert len({r['atom_name'] for r in mapping})==n
    assert 100*287+83*34+83*8==32186
    assert len(manifest['jobs'])==10
    for name in ['PMMA3','PMMA5']:
        jobs=[j for j in manifest['jobs'] if j['molecule']==name]
        assert [j['purpose'] for j in jobs]==['fit','fit','evaluation']


def test_mdp_schedule_and_branch():
    protocol=read_json(ROOT/'config/protocol.json')['stages']
    total=0
    for name,stage in protocol.items():
        text=(ROOT/'mdp'/(name+'.mdp')).read_text()
        parsed={}
        for line in text.splitlines():
            if '=' not in line.split(';')[0]:continue
            k,v=[x.strip() for x in line.split(';')[0].split('=',1)]
            assert k not in parsed;parsed[k]=v
        assert int(parsed['nsteps'])==stage['nsteps']
        assert int(parsed['init-step'])==stage['init_step']
        if not stage['optional']:total+=stage['nsteps']*.001/1000
        if name=='heat':assert parsed['gen-vel']=='yes' and parsed['continuation']=='no'
        else:assert parsed['gen-vel']=='no' and parsed['continuation']=='yes'
        if stage['pressure_bar'] is None:assert parsed['pcoupl']=='no'
        else:assert parsed['pcoupl']=='c-rescale' and float(parsed['ref-p'])==stage['pressure_bar']
        if name=='cool':
            t=np.fromstring(parsed['annealing-time'],sep=' ')
            temp=np.fromstring(parsed['annealing-temp'],sep=' ')
            assert t[0]==stage['start_ps'] and t[-1]==stage['end_ps']
            assert np.allclose(-np.diff(temp)/np.diff(t)*1000,100)
    assert np.isclose(total,77.2)
    assert protocol['nvt_rt']['init_step']==protocol['eq_rt']['init_step']+protocol['eq_rt']['nsteps']
    assert protocol['prod_rt_nvt']['end_ps']==77400
    for path in (ROOT/'mdp').glob('*.mdp'):
        values={}
        for line in path.read_text().splitlines():
            body=line.split(';')[0]
            if '=' in body:
                k,v=[x.strip() for x in body.split('=',1)];values[k]=v
        if 'nstcalcenergy' in values:
            for key in ['nstenergy','nstlog']:
                assert int(values[key])%int(values['nstcalcenergy'])==0


def test_log_rejects_failed_or_unoptimized(tmp_path):
    f=tmp_path/'fixture.log'
    f.write_text('SCF Done: E(RHF) = -1.0\nNormal termination of Gaussian 16\n')
    check_log(f)
    with pytest.raises(ValueError):check_log(f,True)
    f.write_text('SCF Done: x\nError termination\nNormal termination of Gaussian 16\n')
    with pytest.raises(ValueError):check_log(f)


def test_high_precision_gro_round_trip(tmp_path):
    from rdkit import Chem
    from parmed.gromacs import GromacsGroFile
    from forcefield import gro_text
    rows=read_json(ROOT/'structures/OTF_atom_map.json')
    mol=Chem.SDMolSupplier(str(ROOT/'structures/OTF_c01.sdf'),removeHs=False)[0]
    xyz=mol.GetConformer().GetPositions()
    path=tmp_path/'precision.gro';path.write_text(gro_text(rows,xyz,'OTF'))
    loaded=GromacsGroFile.parse(str(path))
    assert [a.name for a in loaded.atoms]==[r['atom_name'] for r in rows]
    assert np.allclose(loaded.coordinates-loaded.coordinates.mean(axis=0),xyz-xyz.mean(axis=0),atol=1e-4)


def test_esp_round_trip_units(tmp_path):
    atoms=np.array([[0.,0.,0.],[2.,0.,0.]])
    points=np.array([[4.,3.,0.],[0.,4.,3.],[-3.,2.,4.]])
    charges=np.array([.3,-.3])
    values=(1/np.linalg.norm(points[:,None]-atoms[None],axis=2))@charges
    path=tmp_path/'synthetic.esp';write_esp(path,atoms,points,values)
    xyz,grid,v=read_esp(path)
    assert np.allclose(xyz,atoms) and np.allclose(grid,points) and np.allclose(v,values)



@pytest.mark.skipif(shutil.which('resp') is None,reason='Native Amber RESP unavailable')
@pytest.mark.parametrize('variant',['OTF_single','OTF_multi','PMMA'])
def test_native_resp_constraints_synthetic_only(tmp_path,variant):
    """Verify actual input parsing/constraints against native RESP, without QM."""
    from rdkit import Chem
    group='PMMA' if variant=='PMMA' else 'OTF'
    if group=='OTF':
        jobs=[dict(id='synthetic_otf_1',molecule='OTF',charge=-1),
              dict(id='synthetic_otf_2',molecule='OTF',charge=-1)]
        if variant=='OTF_single':jobs=jobs[:1]
    else:
        jobs=[dict(id='synthetic_pmma3_1',molecule='PMMA3',charge=0),
              dict(id='synthetic_pmma3_2',molecule='PMMA3',charge=0),
              dict(id='synthetic_pmma5_1',molecule='PMMA5',charge=0),
              dict(id='synthetic_pmma5_2',molecule='PMMA5',charge=0)]
    maps=[maps_for(j) for j in jobs]
    rng=np.random.default_rng(91)
    true={k:float(rng.normal(scale=.2)) for k in sorted({r['equivalence_key'] for m in maps for r in m})}
    if group=='PMMA':
        for template in ['head','internal','tail']:
            rows=maps[0];unit=next(r['unit'] for r in rows if r['template']==template)
            selected=[r for r in rows if r['unit']==unit]
            anchor=next(r for r in selected if r['element']=='C')
            true[anchor['equivalence_key']]-=sum(true[r['equivalence_key']] for r in selected)
    else:
        anchor=next(r for r in maps[0] if r['element']=='C')
        true[anchor['equivalence_key']]+=-1-sum(true[r['equivalence_key']] for r in maps[0])
    esp_paths=[]
    for i,(job,rows) in enumerate(zip(jobs,maps)):
        mol=Chem.SDMolSupplier(str(ROOT/'structures'/(job['molecule']+'_c01.sdf')),removeHs=False)[0]
        atoms=mol.GetConformer().GetPositions()/.529177210903
        # Local, diverse point clouds outside the nuclei. Synthetic Coulomb potential only.
        points=[]
        for center in atoms:
            directions=rng.normal(size=(35,3));directions/=np.linalg.norm(directions,axis=1)[:,None]
            pts=center+directions*5
            pts=pts[np.min(np.linalg.norm(pts[:,None]-atoms[None],axis=2),axis=1)>2]
            points.extend(pts)
        points=np.array(points);q=np.array([true[r['equivalence_key']] for r in rows])
        potential=(1/np.linalg.norm(points[:,None]-atoms[None],axis=2))@q
        path=tmp_path/f'fixture{i}.esp';write_esp(path,atoms,points,potential);esp_paths.append(path)
    (tmp_path/'all.esp').write_text(''.join(p.read_text() for p in esp_paths))
    text,_=controls(group,jobs)
    # Exact-recovery fixture tests format and constraints, not RESP regularization accuracy.
    text=text.replace('qwt=0.0005','qwt=0.0')
    (tmp_path/'fixture.in').write_text(text)
    cp=subprocess.run(['resp','-O','-i','fixture.in','-e','all.esp','-o','fixture.out',
                       '-t','fixture.qout','-p','fixture.pch'],cwd=tmp_path,capture_output=True,text=True)
    out=(tmp_path/'fixture.out').read_text() if (tmp_path/'fixture.out').exists() else cp.stdout+cp.stderr
    if sys.platform=='darwin' and cp.returncode!=0 and 'dyld cache' in out:
        pytest.skip('Host native RESP loader failed before reading input; repeat this native test on Linux.')
    assert cp.returncode==0,out
    # qwt=0 performs a direct ESP solve and does not print an iterative RESP convergence line.
    fitted=read_qout(tmp_path/'fixture.qout',sum(map(len,maps)))
    expected=np.array([true[r['equivalence_key']] for rows in maps for r in rows])
    assert np.max(np.abs(fitted-expected))<3e-4
    if group=='PMMA':
        offset=0
        for rows in maps:
            q=fitted[offset:offset+len(rows)];offset+=len(rows)
            assert abs(q.sum())<1e-5
            for unit in {r['unit'] for r in rows}:
                assert abs(sum(v for r,v in zip(rows,q) if r['unit']==unit))<1e-5
        # Exercise native stage-2 frozen charges and methyl-only equivalence as well.
        text2,_=controls(group,jobs,stage=2)
        (tmp_path/'stage2.in').write_text(text2.replace('qwt=0.001','qwt=0.0'))
        cp2=subprocess.run(['resp','-O','-i','stage2.in','-e','all.esp','-q','fixture.qout',
                            '-o','stage2.out','-t','stage2.qout','-p','stage2.pch'],cwd=tmp_path,capture_output=True,text=True)
        if sys.platform=='darwin' and cp2.returncode!=0 and 'dyld cache' in cp2.stderr:
            pytest.skip('Host native RESP loader failed on stage 2; repeat this native test on Linux.')
        assert cp2.returncode==0,cp2.stdout+cp2.stderr
        q2=read_qout(tmp_path/'stage2.qout',len(fitted))
        mask=np.array([not r['stage2_free'] for rows in maps for r in rows])
        assert np.allclose(q2[mask],fitted[mask],atol=1e-6)
        assert np.max(np.abs(q2-expected))<3e-4


def test_unconfigured_server_rejected_and_render_bash_syntax(tmp_path,monkeypatch):
    import server
    manifest=read_json(ROOT/'structures/manifest.json')
    bundle=tmp_path/'portable_bundle'
    shutil.copytree(ROOT/'structures',bundle/'structures')
    shutil.copytree(ROOT/'gaussian',bundle/'gaussian')
    (bundle/'reports').mkdir()
    monkeypatch.setattr(server,'ROOT',bundle)
    profile=tmp_path/'profile.json'
    cfg=read_json(ROOT/'config/server.example.json')
    profile.write_text(json.dumps(cfg))
    with pytest.raises(ValueError):server.render(profile,'mk')
    assert not (bundle/'runs').exists()
    # Synthetic test allocation, not a statement about available lab resources.
    cfg.update(configured=True,nproc=2,memory_gb=2,pbs_memory_gb=3,queue='TEST_ONLY',
               walltime='01:00:00',scratch_root=str(tmp_path))
    profile.write_text(json.dumps(cfg));server.render(profile,'mk')
    for job in manifest['jobs']:
        folder=bundle/'runs/qm'/job['id']
        assert '@' not in (folder/'opt.gjf').read_text()
        assert 'Pop=MK' in (folder/'esp.gjf').read_text()
        assert 'Geom=AllCheck' in (folder/'esp.gjf').read_text()
        cp=subprocess.run(['bash','-n',str(folder/'job.pbs')],capture_output=True,text=True)
        assert cp.returncode==0,cp.stderr
