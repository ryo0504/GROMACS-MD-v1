"""Prepare/run Amber RESP with explicit multi-conformer/template constraints."""
import argparse
from collections import defaultdict
from pathlib import Path
import re
import numpy as np
from common import ROOT, read_json, write_json, read_esp, esp_error, execute, sha256


def group_jobs(group):
    jobs=read_json(ROOT/'structures/manifest.json')['jobs']
    return [j for j in jobs if j['purpose']=='fit' and
            (j['molecule'] in ['PMMA3','PMMA5'] if group=='PMMA' else j['molecule']==group)]


def maps_for(job, independent=False):
    rows=read_json(ROOT/'structures'/(job['molecule']+'_atom_map.json'))
    if independent:
        # No repeat-template sharing: keep methyl H equivalence but separate each repeat.
        rows=[dict(r,equivalence_key=str(r['unit'])+':'+r['equivalence_key']) for r in rows]
    return rows


def pairs_card(pairs):
    flat=[x for pair in pairs for x in pair]
    return ''.join(''.join(f'{x:5d}' for x in flat[i:i+16])+'\n' for i in range(0,len(flat),16))


def controls(group,jobs,weights=None,stage=1,independent=False):
    """RESP's -99 sentinel omits redundant per-conformer total-charge constraints.

    Actual molecular/template charges are imposed explicitly once on independent
    variable groups; duplicated total constraints would make the KKT system singular.
    The native Amber RESP parser/solver is checked with a small synthetic fixture.
    """
    maps=[maps_for(j,independent) for j in jobs]
    nmol=len(jobs); weights=weights or [1.]*nmol
    shared=defaultdict(list)
    lines=[f'{group} explicit constrained RESP stage {stage}', '', ' &cntrl',
           f' nmol={nmol}, iqopt={1 if stage==1 else 2}, ihfree=1, irstrnt=1,',
           f' qwt={.0005 if stage==1 else .001}, ioutopt=1,', ' &end']
    for im,(job,rows,wt) in enumerate(zip(jobs,maps,weights),1):
        lines += [f'{wt:10.5f}',job['id'],f'{-99:5d}{len(rows):5d}']
        first={}
        for r in rows:
            key=r['equivalence_key']
            frozen=stage==2 and not r['stage2_free']
            if frozen:flag=-1
            elif key in first:flag=first[key]
            else:
                flag=0;first[key]=r['index'];shared[key].append((im,r['index']))
            lines.append(f"{r['atomic_number']:5d}{flag:5d}")
        # For one molecule constraints belong here. For multiple molecules, global section below.
        if nmol==1:
            pairs=[(1,r['index']) for r in rows]
            lines += [f"{len(pairs):5d}{job['charge']:10.5f}",pairs_card(pairs).rstrip('\n')]
        lines.append('')
    constraint_specs=[]
    if nmol>1:
        if group=='PMMA' and not independent:
            for template in ['head','internal','tail']:
                rows=maps[0]
                unit=next(r['unit'] for r in rows if r['template']==template)
                pairs=[(1,r['index']) for r in rows if r['unit']==unit]
                constraint_specs.append((pairs,0.))
        else:
            constraint_specs.append(([(1,r['index']) for r in maps[0]],float(jobs[0]['charge'])))
        for pairs,q in constraint_specs:
            lines += [f'{len(pairs):5d}{q:10.5f}',pairs_card(pairs).rstrip('\n')]
        lines.append('')
        for key,pairs in shared.items():
            if len(pairs)>1:lines += [f'{len(pairs):5d}',pairs_card(pairs).rstrip('\n')]
        lines.append('')
    return '\n'.join(lines)+'\n',maps


def templates():
    for group in ['TPS','OTF','PMMA']:
        jobs=group_jobs(group)
        text,maps=controls(group,jobs)
        folder=ROOT/'resp_templates'/group;folder.mkdir(parents=True,exist_ok=True)
        (folder/'stage1.in.template').write_text(text)
        if group=='PMMA':(folder/'stage2.in.template').write_text(controls(group,jobs,stage=2)[0])
        keys=sorted({r['equivalence_key'] for rows in maps for r in rows})
        constraints=[]
        if group=='PMMA':
            for template in ['head','internal','tail']:
                rows=maps[0];unit=next(r['unit'] for r in rows if r['template']==template)
                row=[sum(r['equivalence_key']==key and r['unit']==unit for r in rows) for key in keys]
                constraints.append(row)
            if np.linalg.matrix_rank(np.array(constraints))!=3:raise ValueError('Template constraints not independent')
        write_json(folder/'constraints.json',dict(group=group,jobs=[j['id'] for j in jobs],
            variable_keys=keys,neutral_template_rows=constraints,
            status='INPUT_TEMPLATE_ONLY; weights regenerated from real ESP point counts',
            note='-99 is an Amber RESP control sentinel, not a physical molecular charge.'))
    print('RESP input templates and constraint maps written. No real charge calculation performed.')


def read_qout(path,n):
    q=np.array([float(x.replace('D','E')) for x in Path(path).read_text().split()])
    if len(q)!=n or not np.isfinite(q).all():raise ValueError('RESP output charge count/values invalid')
    return q


def split_q(q,maps):
    out=[];offset=0
    for rows in maps:
        out.append(q[offset:offset+len(rows)]);offset+=len(rows)
    return out


def run_fit(folder,group,jobs,independent=False):
    folder.mkdir(parents=True,exist_ok=False)
    paths=[ROOT/'runs/qm'/j['id']/'qm.esp' for j in jobs]
    point_counts=[len(read_esp(p)[1]) for p in paths]
    # Native RESP squares wtmol. Equalize point count without renormalizing qwt arbitrarily.
    reference=float(np.mean(point_counts));weights=[np.sqrt(reference/n) for n in point_counts]
    (folder/'combined.esp').write_text(''.join(p.read_text() for p in paths))
    commands=[]
    q=None
    for stage in ([1,2] if group=='PMMA' else [1]):
        text,maps=controls(group,jobs,weights,stage,independent)
        (folder/f'stage{stage}.in').write_text(text)
        args=['resp','-O','-i',f'stage{stage}.in','-e','combined.esp','-o',f'stage{stage}.out',
              '-p',f'stage{stage}.pch','-t',f'stage{stage}.qout','-s',f'stage{stage}.esout']
        if stage==2:args+=['-q','stage1.qout']
        commands.append(execute(args,folder,folder/f'stage{stage}.command.log'))
        result=(folder/f'stage{stage}.out').read_text(errors='replace')
        if not re.search(r'Convergence in\s+\d+ iterations',result):
            raise ValueError('RESP convergence not confirmed; inspect output.')
        if re.search(r'WARNING|FAILED|singular|not enough',result,re.I):
            raise ValueError('RESP output requires review; see warning/failure.')
        q=read_qout(folder/f'stage{stage}.qout',sum(map(len,maps)))
    by_job=split_q(q,maps)
    write_json(folder/'fit_execution.json',dict(commands=commands,point_counts=point_counts,
        native_resp_weights=weights,input_hashes={j['id']:sha256(p) for j,p in zip(jobs,paths)}))
    return maps,by_job


def fit(group):
    jobs=group_jobs(group)
    all_jobs=read_json(ROOT/'structures/manifest.json')['jobs']
    relevant=[j for j in all_jobs if j['molecule'] in ['PMMA3','PMMA5']] if group=='PMMA' else jobs
    methods=set()
    for j in relevant:
        folder=ROOT/'runs/qm'/j['id'];report=read_json(folder/'result.json')
        if report['esp_sha256']!=sha256(folder/'qm.esp'):raise ValueError('ESP changed after extraction')
        if report['optimized_sdf_sha256']!=sha256(folder/'optimized.sdf'):raise ValueError('Optimized geometry changed')
        methods.add(report['esp_method'])
    if len(methods)!=1:raise ValueError('Mixed ESP grids in one fit are prohibited')
    geometry_checks=[]
    if group=='PMMA':
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign
        for name in ['PMMA3','PMMA5']:
            selected=[j for j in relevant if j['molecule']==name]
            for i,j in enumerate(selected):
                for k in selected[:i]:
                    a=Chem.RemoveHs(Chem.SDMolSupplier(str(ROOT/'runs/qm'/j['id']/'optimized.sdf'),removeHs=False)[0])
                    b=Chem.RemoveHs(Chem.SDMolSupplier(str(ROOT/'runs/qm'/k['id']/'optimized.sdf'),removeHs=False)[0])
                    rms=float(rdMolAlign.GetBestRMS(a,b,maxMatches=100))
                    geometry_checks.append(dict(first=j['id'],second=k['id'],heavy_atom_rmsd_A=rms))
                    if rms<.1:raise ValueError('Optimized PMMA conformers collapsed to a duplicate (<0.1 Å). Generate a new independent conformer before fitting/evaluation.')
    folder=ROOT/'runs/resp'/group
    maps,qs=run_fit(folder,group,jobs)
    templates_q={}
    for rows,q in zip(maps,qs):
        for r,value in zip(rows,q):
            key=r['equivalence_key']
            if key in templates_q and abs(templates_q[key]-value)>1e-6:raise ValueError('RESP equivalence mismatch')
            templates_q[key]=float(value)
    validation=[]
    for j in relevant:
        rows=maps_for(j);q=np.array([templates_q[r['equivalence_key']] for r in rows])
        if abs(q.sum()-j['charge'])>1e-5:raise ValueError(f"Net charge mismatch: {j['id']}; inspect RESP rounding/fit")
        validation.append(dict(job=j['id'],purpose=j['purpose'],
            **esp_error(ROOT/'runs/qm'/j['id']/'qm.esp',q)))
    targets=['PMMA19'] if group=='PMMA' else [group]
    charges_folder=ROOT/'runs/resp/charges';charges_folder.mkdir(parents=True,exist_ok=True)
    if group=='PMMA':
        # Independent molecular RESP baseline, trained on the same fit conformers.
        for molecule in ['PMMA3','PMMA5']:
            subset=[j for j in jobs if j['molecule']==molecule]
            independent_maps,independent_q=run_fit(folder/('independent_'+molecule),'PMMA',subset,independent=True)
            baseline={}
            for r,value in zip(independent_maps[0],independent_q[0]):baseline[r['equivalence_key']]=float(value)
            for entry in validation:
                j=next(x for x in relevant if x['id']==entry['job'])
                if j['molecule']!=molecule:continue
                q=np.array([baseline[r['equivalence_key']] for r in maps_for(j,True)])
                base=esp_error(ROOT/'runs/qm'/j['id']/'qm.esp',q)
                entry['independent_rmse_au']=base['rmse_au']
                entry['rmse_ratio']=entry['rmse_au']/max(base['rmse_au'],1e-12)
                entry['transfer_review_required']=entry['rmse_ratio']>1.25
    final_hashes={}
    for molecule in targets:
        rows=read_json(ROOT/'structures'/(molecule+'_atom_map.json'))
        q=np.array([templates_q[r['equivalence_key']] for r in rows])
        charge=0 if molecule.startswith('PMMA') else (1 if molecule=='TPS' else -1)
        if abs(q.sum()-charge)>1e-5:
            raise ValueError('Transferred charge sum exceeds 1e-5 e. Do not apply an arbitrary offset; review native RESP precision.')
        path=charges_folder/(molecule+'.charges')
        if path.exists():raise ValueError('Charge file exists; preserve the adopted input')
        # More printed digits do not add precision to native RESP qout.
        np.savetxt(path,q,fmt='%.10f');final_hashes[molecule]=sha256(path)
    write_json(folder/'review.json',dict(status='CALCULATED_REVIEW_REQUIRED',
        approved=False,esp_method=next(iter(methods)),charges_sha256=final_hashes,
        template_charges=templates_q,validation=validation,geometry_checks=geometry_checks,
        checks_to_review=['charge/geometry/order','ESP scatter and residuals','held-out conformers',
                          'template transfer and local neutrality','rounding precision','missing force-field parameters']))
    print('Real RESP completed. Inspect',folder/'review.json','before adopting the force field.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['templates','fit']);p.add_argument('--group',choices=['TPS','OTF','PMMA'])
    a=p.parse_args()
    if a.action=='templates':templates()
    elif not a.group:p.error('--group is required for fit')
    else:fit(a.group)
