from pathlib import Path
import json
import numpy as np
import pytest
from mdworkbench.analysis.ethanol_rdf import production_indices, intermolecular_rdf
from mdworkbench.ethanol_benchmark import config, build, scientific_inputs, mdp

CONFIG = Path(__file__).resolve().parents[1] / 'systems/ethanol_27/system.yaml'


def test_frames():
    steps = np.arange(0, 601, 5)
    take = production_indices(steps)
    assert len(take) == 80
    assert steps[take[0]] == 205 and steps[take[-1]] == 600
    assert 200 not in steps[take]
    for invalid in (steps[:-1], np.append(steps,600), steps[::-1]):
        with pytest.raises(ValueError): production_indices(invalid)


def test_pair_exclusion_mic_normalization():
    # 27 centers at x=.1; 27 H at x=13.7. MIC O-H distance=.2 A.
    xyz = np.zeros((1,54,3)); xyz[0,:27,0]=.1; xyz[0,27:,0]=13.7
    ids = np.tile(np.arange(27),2)
    result = intermolecular_rdf(xyz,np.array([[13.8]*3]),ids,np.arange(27),np.arange(27,54))
    assert len(result['OO']) == 45
    assert result['OH_counts'].sum() == 27*26
    assert result['OO_counts'][0] == 27*26
    index=np.searchsorted(result['edges_A'],.2)-1
    shell=4*np.pi/3*np.diff(result['edges_A']**3)
    expected=(27*26)/(27*26/13.8**3*shell[index])
    assert result['OH'][index] == pytest.approx(expected)
    # Independent distance histogram oracle, including a pair beyond r_max.
    xyz[0,27:,0]=7.0
    result=intermolecular_rdf(xyz,np.array([[13.8]*3]),ids,np.arange(27),np.arange(27,54))
    assert result['OH_counts'].sum() == 0


def test_no_unreviewed_science():
    cfg=config(CONFIG)
    with pytest.raises(ValueError,match='GAFF'): scientific_inputs(cfg,CONFIG.parent)
    cfg['md']['steps']=601
    import yaml, tempfile
    with tempfile.TemporaryDirectory() as temporary:
        path=Path(temporary)/'bad.yaml'; path.write_text(yaml.safe_dump(cfg))
        with pytest.raises(ValueError,match='steps'): config(path)


def test_build_reproducible_velocity_and_density(tmp_path):
    pytest.importorskip('ase'); pytest.importorskip('rdkit')
    cfg=config(CONFIG)
    a=build(cfg,tmp_path/'a'); b=build(cfg,tmp_path/'b')
    x=np.load(a/'initial.npz'); y=np.load(b/'initial.npz')
    assert x['positions_A'].shape==(243,3)
    assert np.array_equal(x['positions_A'],y['positions_A'])
    assert np.array_equal(x['velocities_nm_ps'],y['velocities_nm_ps'])
    report=json.loads((a/'build_report.json').read_text())
    assert report['temperature_726_dof_K']==pytest.approx(300)
    assert report['temperature_729_dof_K']==pytest.approx(300*726/729)
    assert report['density_g_cm3']==pytest.approx(.786,abs=.001)
    assert np.linalg.norm(report['com_momentum_ase'])<1e-12
    # Molecule COM velocities remain nonzero.
    momentum=x['velocities_nm_ps']*x['masses_amu'][:,None]
    assert np.linalg.norm(momentum.reshape(27,9,3).sum(axis=1))>0


def test_mdp_not_generic_staged_protocol():
    cfg=config(CONFIG)
    cfg['nonbonded'].update(coulombtype='PME',vdwtype='Cut-off',rlist_nm=.6,rcoulomb_nm=.6,rvdw_nm=.6,dispersion_correction='no')
    values=dict(line.strip().split(' = ') for line in mdp(cfg).splitlines())
    assert values['integrator']=='md-vv' and values['nsteps']=='600'
    assert values['gen-vel']=='no' and values['pcoupl']=='no'
    assert values['constraints']=='none' and values['nh-chain-length']=='3'
    assert values['comm-mode']=='none' and values['nstxout-compressed']=='5'


def test_prepare_preserves_coordinates_velocities_and_rejects_tampering(tmp_path, monkeypatch):
    """Synthetic zero-charge topology tests file mechanics, NOT a RESP/GAFF result."""
    pmd=pytest.importorskip('parmed')
    pytest.importorskip('rdkit')
    mda=pytest.importorskip('MDAnalysis')
    from ase.data import atomic_numbers
    from mdworkbench import ethanol_benchmark as bench
    import yaml
    monkeypatch.setattr(bench,'collect_metadata',lambda *args,**kwargs: {'commands':[]})
    cfg=config(CONFIG)
    initial=bench.build(cfg,tmp_path/'initial')
    report=json.loads((initial/'build_report.json').read_text())
    data=np.load(initial/'initial.npz')
    single=pmd.Structure()
    for i,symbol in enumerate(report['ase_single_symbols']):
        atom=pmd.Atom(name=f'{symbol}{i}',type=symbol,mass=data['masses_amu'][i],charge=0,atomic_number=atomic_numbers[symbol])
        atom.atom_type=pmd.AtomType(symbol,i+1,atom.mass,atomic_number=atom.atomic_number)
        atom.atom_type.set_lj_params(.1,1.5)
        single.add_atom(atom,'ETH',1)
    bond_type=pmd.BondType(300,1.1)
    single.bond_types.append(bond_type)
    for i,j in report['ase_single_bonds']:
        single.bonds.append(pmd.Bond(single.atoms[i],single.atoms[j],type=bond_type))
    single.save(str(tmp_path/'synthetic.top'))
    # Supply every graph-distance-3 pair; no implicit missing-1-4 fallback.
    graph={i:set() for i in range(9)}
    for i,j in report['ase_single_bonds']:
        graph[i].add(j); graph[j].add(i)
    pairs=[]
    for i in range(9):
        visited={i}; frontier={i}
        for depth in range(1,4):
            frontier=set().union(*(graph[j] for j in frontier))-visited
            visited |= frontier
        pairs.extend((i,j) for j in frontier if i<j)
    fixture=tmp_path/'synthetic.top'
    text=fixture.read_text()
    text=text.replace('[ system ]','[ pairs ]\n'+'\n'.join(f'{i+1} {j+1} 1' for i,j in pairs)+'\n\n[ system ]')
    fixture.write_text(text)
    (tmp_path/'synthetic.charges').write_text('0 '*9)
    cfg['force_field']['family']='gaff2'
    cfg['resp'].update(qm_backend='SYNTHETIC_TEST_ONLY',method='TEST',basis='TEST',geometry_file='initial/ethanol_single.xyz',
        esp_conditions='SYNTHETIC_TEST_ONLY',fitting_conditions='NOT_REAL_RESP',source_file='synthetic.charges',atom_order_confirmed=True)
    cfg['structure']['ase_atom_order']=list(range(9))
    cfg['nonbonded'].update(coulombtype='PME',vdwtype='Cut-off',rlist_nm=.6,rcoulomb_nm=.6,rvdw_nm=.6,
        dispersion_correction='no',adoption_reason='TEST_FIXTURE_ONLY')
    cfg['review'].update(force_field_parameters_confirmed=True,max_force_kJ_mol_nm=10000)
    path=tmp_path/'system.yaml'; path.write_text(yaml.safe_dump(cfg))
    prepared=bench.prepare(cfg,path,initial,tmp_path/'synthetic.top',tmp_path/'run')
    loaded=pmd.load_file(str(prepared/'topology/system.top'))
    assert len(loaded.atoms)==243 and len(loaded.residues)==27
    u=mda.Universe(str(prepared/'initial.gro'))
    assert np.allclose(u.atoms.positions,data['positions_A'],atol=1e-5)
    # MDAnalysis default velocities are A/ps.
    assert np.allclose(u.trajectory.ts.velocities/10,data['velocities_nm_ps'],atol=1e-6)
    bench.verify_snapshot(prepared)
    (prepared/'nvt.mdp').write_text('changed')
    with pytest.raises(ValueError,match='changed'): bench.verify_snapshot(prepared)


def test_thermodynamics_uses_labels_not_selection_order(tmp_path):
    from mdworkbench.ethanol_benchmark import read_thermodynamics
    path=tmp_path/'energy.xvg'
    path.write_text('@ s0 legend "Potential"\n@ s1 legend "Temperature"\n0 -123 298.765\n0.3 -124 301\n')
    time,temperature,potential=read_thermodynamics(path)
    assert np.array_equal(temperature,[298.765,301])
    assert np.array_equal(potential,[-123,-124])
    assert np.array_equal(time,[0,.3])
