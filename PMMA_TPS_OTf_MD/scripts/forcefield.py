"""Build actual GAFF2 inputs only after reviewed real RESP output exists."""
import argparse
from pathlib import Path
import re
import shutil
import numpy as np
from common import ROOT, read_json, write_json, sha256, execute


def sections(text):
    blocks=[]; current=None;lines=[]
    for line in text.splitlines():
        m=re.match(r'\s*\[\s*([\w]+)\s*\]',line)
        if m:
            if current is not None:blocks.append((current,lines))
            current=m.group(1);lines=[line]
        elif current is not None:lines.append(line)
    if current is not None:blocks.append((current,lines))
    return blocks


def body(lines):
    return [l.split(';',1)[0].strip() for l in lines[1:] if l.split(';',1)[0].strip()]


def gro_text(rows,xyz_A,resname,precision=5):
    xyz_nm=np.asarray(xyz_A)/10.;xyz_nm=xyz_nm-xyz_nm.mean(axis=0)+15
    lines=['GAFF2 RESP; coordinates centered in 30 nm validation box',str(len(rows))]
    width=precision+5
    for r,(x,y,z) in zip(rows,xyz_nm):
        resid=r.get('resid',1);rn=r.get('resname',resname)
        lines.append(f"{resid:5d}{rn:<5}{r['atom_name']:>5}{r['index']:5d}"
                     f'{x:{width}.{precision}f}{y:{width}.{precision}f}{z:{width}.{precision}f}')
    return '\n'.join(lines)+'\n  30.00000  30.00000  30.00000\n'


def build(name):
    import parmed as pmd
    from rdkit import Chem
    group='PMMA' if name=='PMMA19' else name
    review=read_json(ROOT/'runs/resp'/group/'review.json')
    if review.get('approved') is not True:raise ValueError('Review real ESP/RESP/transfer results, then mark review.json approved=true.')
    charge_path=ROOT/'runs/resp/charges'/(name+'.charges')
    if review['charges_sha256'][name]!=sha256(charge_path):raise ValueError('RESP charges changed after review')
    q=np.loadtxt(charge_path,ndmin=1)
    rows=read_json(ROOT/'structures'/(name+'_atom_map.json'))
    charge={'PMMA19':0,'TPS':1,'OTF':-1}[name]
    if len(q)!=len(rows) or abs(q.sum()-charge)>1e-5:raise ValueError('RESP charge number/sum invalid')
    source=ROOT/'structures/PMMA19_c01.sdf' if name=='PMMA19' else ROOT/'runs/qm'/(name+'_c01')/'optimized.sdf'
    mol=Chem.SDMolSupplier(str(source),removeHs=False)[0]
    folder=ROOT/'runs/forcefield'/name;folder.mkdir(parents=True,exist_ok=False)
    shutil.copy2(source,folder/'input.sdf');shutil.copy2(charge_path,folder/'resp.charges')
    commands=[]
    commands.append(execute(['antechamber','-i','input.sdf','-fi','sdf','-o','charged.mol2','-fo','mol2',
        '-at','gaff2','-c','rc','-cf','resp.charges','-nc',str(charge),'-m','1','-s','2'],folder,folder/'antechamber.command.log'))
    amber_mol=pmd.load_file(str(folder/'charged.mol2'))
    if len(amber_mol.atoms)!=len(rows):raise ValueError('Antechamber atom count changed')
    if not np.allclose(amber_mol.coordinates,mol.GetConformer().GetPositions(),atol=2e-4):
        raise ValueError('Antechamber coordinate/order changed; mapping review required')
    if not np.allclose([a.charge for a in amber_mol.atoms],q,atol=2e-6):raise ValueError('Antechamber changed/reordered charges')
    expected={tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    actual={tuple(sorted((b.atom1.idx,b.atom2.idx))) for b in amber_mol.bonds}
    if expected!=actual:raise ValueError('Antechamber changed the bond graph')
    commands.append(execute(['parmchk2','-i','charged.mol2','-f','mol2','-o','parameters.frcmod','-s','gaff2'],folder,folder/'parmchk2.command.log'))
    frcmod=(folder/'parameters.frcmod').read_text()
    if re.search(r'ATTN|NEED REVISION',frcmod,re.I):raise ValueError('Unreviewed force-field parameters in frcmod; stop and inspect.')
    (folder/'tleap.in').write_text('source leaprc.gaff2\nloadamberparams parameters.frcmod\n'
        'MOL = loadmol2 charged.mol2\ncheck MOL\nsaveamberparm MOL system.prmtop system.inpcrd\nquit\n')
    commands.append(execute(['tleap','-f','tleap.in'],folder,folder/'tleap.command.log'))
    amber=pmd.load_file(str(folder/'system.prmtop'),xyz=str(folder/'system.inpcrd'))
    if len(amber.atoms)!=len(rows) or abs(sum(a.charge for a in amber.atoms)-charge)>1e-5:
        raise ValueError('Amber topology charge/count mismatch')
    if any(x.type is None for x in list(amber.bonds)+list(amber.angles)+list(amber.dihedrals)):
        raise ValueError('Missing Amber bonded parameters')
    commands.append(execute(['acpype','-p','system.prmtop','-x','system.inpcrd','-b',name],folder,folder/'acpype.command.log'))
    paths=list(folder.rglob('*_GMX.top'))
    if len(paths)!=1:raise ValueError('Expected one ACPYPE GROMACS topology')
    blocks=sections(paths[0].read_text());resname='PMM' if name=='PMMA19' else name
    output=[]
    for kind,lines in blocks:
        if kind=='moleculetype':
            rows_body=body(lines)
            if len(rows_body)!=1:raise ValueError('Unexpected molecule definitions')
            lines=['[ moleculetype ]',f'{name} {rows_body[0].split()[1]}']
        elif kind=='atoms':
            parsed=body(lines)
            if len(parsed)!=len(rows):raise ValueError('GROMACS atom count mismatch')
            lines=['[ atoms ]','; nr type resnr residue atom cgnr charge mass']
            for fields,r,value in zip([x.split() for x in parsed],rows,q):
                if int(fields[0])!=r['index']:raise ValueError('GROMACS order mismatch')
                fields[3]=resname;fields[4]=r['atom_name'];fields[6]=f'{value:.10f}'
                lines.append(' '.join(fields))
        elif kind=='system':lines=['[ system ]',name+' isolated component-energy validation']
        elif kind=='molecules':lines=['[ molecules ]',name+' 1']
        output.append((kind,lines))
    (folder/'isolated.top').write_text('\n\n'.join('\n'.join(v) for k,v in output)+'\n')
    (folder/(name+'.gro')).write_text(gro_text(rows,amber.coordinates,resname))
    # Match both energy engines to the precision actually written to GRO.
    from parmed.gromacs import GromacsGroFile
    same=GromacsGroFile.parse(str(folder/(name+'.gro')))
    amber.coordinates=same.coordinates
    amber.save(str(folder/'same_coordinates.inpcrd'),overwrite=False)
    write_json(folder/'build_report.json',dict(status='BUILT_ENERGY_AND_PARAMETER_REVIEW_PENDING',
        molecule=name,commands=commands,resp_sha256=sha256(charge_path),
        prmtop_sha256=sha256(folder/'system.prmtop'),top_sha256=sha256(folder/'isolated.top'),
        gro_sha256=sha256(folder/(name+'.gro')),parameter_review_approved=False,
        counts=dict(atoms=len(amber.atoms),bonds=len(amber.bonds),angles=len(amber.angles),dihedrals=len(amber.dihedrals)),
        charge_e=sum(a.charge for a in amber.atoms)))
    print('Generated real',name,'GAFF2 topology. Run component-energy checks and inspect parameters.')


def assemble():
    destination=ROOT/'runs/adopted_input'
    if destination.exists():raise ValueError('Adopted input exists; never overwrite a scientific input set')
    all_types={};defaults=None;molecules={};inputs={};composition={}
    for name in ['PMMA19','TPS','OTF']:
        folder=ROOT/'runs/forcefield'/name
        report=read_json(folder/'build_report.json');energy=read_json(folder/'energy_comparison.json')
        if report.get('parameter_review_approved') is not True or energy.get('passed') is not True:
            raise ValueError(f'Parameter/energy checks are incomplete: {name}')
        for filename,key in [('system.prmtop','prmtop_sha256'),('isolated.top','top_sha256'),(name+'.gro','gro_sha256')]:
            if sha256(folder/filename)!=report[key]:raise ValueError('Topology changed after checks')
        if energy['top_sha256']!=sha256(folder/'isolated.top'):raise ValueError('Energy check for a different topology')
        blocks=sections((folder/'isolated.top').read_text());local=[];inside=False
        for kind,lines in blocks:
            if kind=='defaults':
                values=body(lines)
                if defaults is not None and defaults!=values:raise ValueError('Incompatible global force-field defaults')
                defaults=values
            elif kind=='atomtypes':
                for line in body(lines):
                    fields=line.split();key=fields[0]
                    if key in all_types and all_types[key]!=fields:raise ValueError('Conflicting atom type: '+key)
                    all_types[key]=fields
            elif kind=='moleculetype':inside=True;local+=lines
            elif kind in ['system','molecules']:inside=False
            elif inside:local+=lines
            elif kind not in ['defaults','atomtypes']:raise ValueError('Unrecognized global topology section: '+kind)
        molecules[name]='\n'.join(local)+'\n';inputs[name]=report
        atom_rows=[line.split() for kind,lines in sections(molecules[name]) if kind=='atoms' for line in body(lines)]
        composition[name]=dict(atoms=len(atom_rows),charge_e=sum(float(r[6]) for r in atom_rows),
                               mass_amu=sum(float(r[7]) for r in atom_rows))
    counts=read_json(ROOT/'config/system.json')['counts']
    atoms=sum(counts[name]*c['atoms'] for name,c in composition.items())
    net_charge=sum(counts[name]*c['charge_e'] for name,c in composition.items())
    if atoms!=32186 or abs(net_charge)>1e-3:
        raise ValueError(f'Mixed topology count/charge mismatch: {atoms} atoms, {net_charge} e')
    destination.mkdir(parents=True)
    (destination/'forcefield.itp').write_text('[ defaults ]\n'+'\n'.join(defaults)+'\n')
    (destination/'atomtypes.itp').write_text('[ atomtypes ]\n'+'\n'.join(' '.join(v) for v in all_types.values())+'\n')
    for name,text in molecules.items():
        (destination/(name+'.itp')).write_text(text)
        shutil.copy2(ROOT/'runs/forcefield'/name/(name+'.gro'),destination/(name+'.gro'))
    shutil.copy2(ROOT/'topology/topol.top.in',destination/'topol.top')
    # Strip the template-only comment after all included real files have been generated.
    top=destination/'topol.top';top.write_text('\n'.join(top.read_text().splitlines()[1:])+'\n')
    write_json(destination/'adoption.json',dict(status='COMPONENTS_ADOPTED_MIXED_SYSTEM_PREFLIGHT_PENDING',
        input_reports=inputs,composition=composition,total_atoms=atoms,total_charge_e=net_charge,
        total_mass_amu=sum(counts[name]*c['mass_amu'] for name,c in composition.items()),
        file_hashes={p.name:sha256(p) for p in destination.iterdir() if p.is_file()}))
    print('Adopted component input:',destination,'; mixed-system packing/preflight still required.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('build');a.add_argument('molecule',choices=['PMMA19','TPS','OTF'])
    s.add_parser('assemble');a=p.parse_args()
    if a.action=='build':build(a.molecule)
    else:assemble()
