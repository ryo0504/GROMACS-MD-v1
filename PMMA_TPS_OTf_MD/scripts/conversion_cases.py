"""Prepare a full component-energy matrix from actual GAFF2/RESP topologies."""
import argparse
from pathlib import Path
import shutil
import numpy as np
from common import ROOT, read_json, write_json, sha256
from forcefield import sections,body,gro_text


def combined_top(names):
    defaults=None;types={};molecules=[]
    for name in names:
        blocks=sections((ROOT/'runs/forcefield'/name/'isolated.top').read_text())
        inside=False;local=[]
        for kind,lines in blocks:
            if kind=='defaults':
                d=body(lines)
                if defaults is not None and defaults!=d:raise ValueError('Default mismatch')
                defaults=d
            elif kind=='atomtypes':
                for line in body(lines):
                    f=line.split()
                    if f[0] in types and types[f[0]]!=f:raise ValueError('Atom type conflict')
                    types[f[0]]=f
            elif kind=='moleculetype':inside=True;local+=lines
            elif kind in ['system','molecules']:inside=False
            elif inside:local+=lines
        molecules.append('\n'.join(local))
    return ('[ defaults ]\n'+'\n'.join(defaults)+'\n\n[ atomtypes ]\n'+
            '\n'.join(' '.join(v) for v in types.values())+'\n\n'+
            '\n\n'.join(molecules)+'\n\n[ system ]\nComponent energy test case\n\n[ molecules ]\n'+
            '\n'.join(n+' 1' for n in names)+'\n')


def prepare():
    import parmed as pmd
    from parmed.gromacs import GromacsGroFile
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms
    root=ROOT/'runs/conversion_cases';root.mkdir(parents=True,exist_ok=False)
    structures={};maps={}
    for name in ['PMMA19','TPS','OTF']:
        folder=ROOT/'runs/forcefield'/name
        structures[name]=pmd.load_file(str(folder/'system.prmtop'),xyz=str(folder/'system.inpcrd'))
        maps[name]=read_json(ROOT/'structures'/(name+'_atom_map.json'))
    cases=[]
    def save(name,names,structure,coords,rows,description):
        folder=root/name;folder.mkdir()
        structure.coordinates=coords
        (folder/'isolated.top').write_text(combined_top(names))
        (folder/'coordinates.gro').write_text(gro_text(rows,coords,'TST'))
        structure.coordinates=GromacsGroFile.parse(str(folder/'coordinates.gro')).coordinates
        structure.save(str(folder/'system.prmtop'),overwrite=False)
        structure.save(str(folder/'same_coordinates.inpcrd'),overwrite=False)
        write_json(folder/'case.json',dict(name=name,description=description,
            physical_production_initial_state=False,
            hashes={p.name:sha256(p) for p in folder.iterdir() if p.is_file()}))
        cases.append(str(folder.relative_to(ROOT)))
    mol=Chem.SDMolSupplier(str(ROOT/'structures/PMMA19_c01.sdf'),removeHs=False)[0]
    def atom(unit,role):return next(r['index']-1 for r in maps['PMMA19'] if r['unit']==unit and r['role']==role)
    indices=[atom(4,'backbone_substituted_C'),atom(5,'backbone_head_C'),
             atom(5,'backbone_substituted_C'),atom(6,'backbone_head_C')]
    for delta in [-60,60]:
        m=Chem.Mol(mol)
        conf=m.GetConformer()
        current=rdMolTransforms.GetDihedralDeg(conf,*indices)
        rdMolTransforms.SetDihedralDeg(conf,*indices,current+delta)
        save('PMMA19_torsion_'+('minus60' if delta<0 else 'plus60'),['PMMA19'],
             structures['PMMA19'].copy(pmd.Structure),conf.GetPositions(),maps['PMMA19'],
             'Deterministic backbone torsion change for conversion validation; not an equilibrated/minimized state')
    tps=structures['TPS'];otf=structures['OTF']
    st=next(r['index']-1 for r in maps['TPS'] if r['element']=='S')
    so=next(r['index']-1 for r in maps['OTF'] if r['element']=='S')
    for distance_nm in [.6,.8,1.0]:
        for degrees in [0,90]:
            theta=np.deg2rad(degrees)
            rotation=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
            x=tps.coordinates-tps.coordinates[st]
            y=(otf.coordinates-otf.coordinates[so])@rotation.T+np.array([distance_nm*10,0,0])
            combined=tps+otf
            rows=[]
            for resid,name in enumerate(['TPS','OTF'],1):
                for r in maps[name]:rows.append(dict(r,index=len(rows)+1,resid=resid,resname=name))
            save(f'pair_{distance_nm:.1f}nm_{degrees}deg',['TPS','OTF'],combined,np.concatenate([x,y]),rows,
                 'TPS S--OTf S separation/orientation scan for energy conversion; not an equilibrium ion-pair prediction')
    write_json(root/'matrix.json',dict(cases=cases,status='PREPARED_NOT_CALCULATED',
        additional_required_single_molecules=['PMMA19','TPS','OTF']))
    print('Prepared',len(cases),'additional energy cases. No energies calculated.')


def summarize():
    matrix=read_json(ROOT/'runs/conversion_cases/matrix.json')
    passed=True;evidence=[]
    for name in matrix['additional_required_single_molecules']:
        folder=ROOT/'runs/forcefield'/name
        build=read_json(folder/'build_report.json')
        if build.get('parameter_review_approved') is not True:passed=False
        folders=[folder] if name=='PMMA19' else [folder]
        for f in folders:
            report=read_json(f/'energy_comparison.json')
            passed=passed and report['passed']
            if sha256(f/'isolated.top')!=report['top_sha256']:raise ValueError('Changed tested topology')
            gro=name+'.gro'
            if sha256(f/gro)!=report['gro_sha256']:raise ValueError('Changed tested geometry')
            if sha256(f/'system.prmtop')!=report['prmtop_sha256']:raise ValueError('Changed Amber topology')
            evidence.append(dict(path=str((f/'energy_comparison.json').relative_to(ROOT)),sha256=sha256(f/'energy_comparison.json')))
    for path in matrix['cases']:
        f=ROOT/path;report=read_json(f/'energy_comparison.json')
        passed=passed and report['passed']
        for filename,key in [('isolated.top','top_sha256'),('coordinates.gro','gro_sha256'),('system.prmtop','prmtop_sha256')]:
            if sha256(f/filename)!=report[key]:raise ValueError('Energy case changed')
        evidence.append(dict(path=str((f/'energy_comparison.json').relative_to(ROOT)),sha256=sha256(f/'energy_comparison.json')))
    write_json(ROOT/'reports/conversion_review.json',dict(approved=bool(passed),single_molecules=bool(passed),
        multiple_polymer_geometries=bool(passed),ion_pair_distances_and_orientations=bool(passed),
        status='PASSED' if passed else 'FAILED_OR_PARAMETER_REVIEW_PENDING',evidence_files=evidence))
    if not passed:raise ValueError('Conversion/parameter review matrix incomplete')
    print('All planned conversion tests and parameter reviews passed; mixed-system MD preflight remains.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','summarize'])
    a=p.parse_args();prepare() if a.action=='prepare' else summarize()
