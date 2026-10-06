"""Dedicated no-minimization, continuous 600-step ethanol benchmark."""
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import sys
import warnings
import numpy as np
import yaml
from mdworkbench.provenance.metadata import collect_metadata, hash_file, write_metadata
from mdworkbench.gromacs.prepare import snapshot_topology
from mdworkbench.analysis.ethanol_rdf import production_indices, intermolecular_rdf, read_gromacs, read_mace


def config(path):
    cfg = yaml.safe_load(Path(path).read_text())
    expected = {'system': {'molecule_count': 27, 'atoms_per_molecule': 9},
                'box': {'length_A': 13.8, 'pbc': [True, True, True]},
                'placement': {'grid': [3, 3, 3], 'spacing_A': 4.6, 'rotation_deg': [-5., 5.], 'rng_seed': 42},
                'velocities': {'temperature_K': 300., 'rng_seed': 42, 'degrees_of_freedom': 726},
                'md': {'ensemble': 'NVT', 'temperature_K': 300., 'timestep_fs': .5, 'steps': 600},
                'thermostat': {'chain_length': 3, 'ase_tdamp_fs': 50., 'gromacs_tau_t_ps': .3141592654},
                'trajectory': {'stride_steps': 5},
                'analysis': {'production_start_exclusive_step': 200, 'r_max_A': 6.762, 'rdf_bins': 45}}
    for section, values in expected.items():
        for key, value in values.items():
            if cfg.get(section, {}).get(key) != value:
                raise ValueError(f'Benchmark condition changed: {section}.{key} must be {value}')
    return cfg


def build(cfg, output):
    from ase.build import molecule
    from ase.io import write
    from ase import units
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    single = molecule('CH3CH2OH')
    single.translate(-single.get_center_of_mass())
    single.rotate('x', (1, 1, 1), center=(0, 0, 0))
    rng = np.random.default_rng(42)
    atoms = single.copy(); del atoms[:]
    rotations = []
    for index in np.ndindex(3, 3, 3):
        copy = single.copy()
        angles = rng.uniform(-5, 5, 3)
        rotations.append(angles.tolist())
        for axis, angle in zip('xyz', angles):
            copy.rotate(angle, axis, center=(0, 0, 0))
        copy.translate((np.asarray(index) + .5) * 4.6)
        atoms += copy
    atoms.set_cell([13.8] * 3); atoms.set_pbc(True)
    mass = atoms.get_masses()
    # Independent RNG: velocity seed 42 is reset after the placement RNG.
    p = np.random.default_rng(42).standard_normal((243, 3)) * np.sqrt(mass[:, None] * units.kB * 300)
    p -= mass[:, None] * (p.sum(axis=0) / mass.sum())
    kinetic = .5 * np.sum(p * p / mass[:, None])
    p *= np.sqrt(726 * units.kB * 300 / (2 * kinetic))
    atoms.set_momenta(p)
    write(root / 'ethanol_single.xyz', single)
    write(root / 'initial.extxyz', atoms)
    np.savez(root / 'initial.npz', positions_A=atoms.positions, velocities_nm_ps=atoms.get_velocities() * units.fs * 100,
             masses_amu=mass, molecule_ids=np.repeat(np.arange(27), 9))
    # Explicit bond graph and ASE atom order, without RDKit geometry optimization.
    from rdkit import Chem
    rw = Chem.RWMol()
    for symbol in single.get_chemical_symbols():
        rw.AddAtom(Chem.Atom(symbol))
    bonds = []
    for i in range(9):
        for j in range(i + 1, 9):
            d = np.linalg.norm(single.positions[i] - single.positions[j])
            if ((single[i].symbol == 'H') != (single[j].symbol == 'H') and d < 1.25) or (single[i].symbol != 'H' and single[j].symbol != 'H' and d < 1.65):
                rw.AddBond(i, j, Chem.BondType.SINGLE); bonds.append((i, j))
    mol = rw.GetMol(); Chem.SanitizeMol(mol)
    conformer = Chem.Conformer(9)
    for i, xyz in enumerate(single.positions): conformer.SetAtomPosition(i, xyz.tolist())
    mol.AddConformer(conformer)
    with Chem.SDWriter(str(root / 'ethanol_single.sdf')) as writer: writer.write(mol)
    oxygen = single.get_chemical_symbols().index('O')
    h = [j if i == oxygen else i for i, j in bonds if oxygen in (i,j) and single[j if i == oxygen else i].symbol == 'H']
    if len(bonds) != 8 or len(h) != 1: raise ValueError('Invalid ethanol bond graph')
    density = mass.sum() * 1.66053906660 / atoms.get_volume()
    report = {'atom_count': 243, 'molecule_count': 27, 'molecular_mass_amu': single.get_masses().sum(),
              'box_A': [13.8]*3, 'volume_A3': atoms.get_volume(), 'density_g_cm3': density,
              'density_warning': bool(abs(density - .786) > .01), 'rotations_deg': rotations,
              'ase_single_symbols': single.get_chemical_symbols(), 'ase_single_bonds': bonds,
              'oxygen_ase_index': oxygen, 'hydroxyl_h_ase_index': h[0],
              'com_momentum_ase': atoms.get_momenta().sum(axis=0).tolist(),
              'temperature_726_dof_K': 2*atoms.get_kinetic_energy()/(726*units.kB),
              'temperature_729_dof_K': atoms.get_temperature(),
              'physical_com_removed_dof': 726, 'ase_unconstrained_dof': 729,
              'gromacs_comm_none_dof': 729, 'reference_initial_state_verified': False}
    write_metadata(root / 'build_report.json', report)
    metadata = collect_metadata([root/'initial.npz', root/'initial.extxyz'], 42, cfg)
    metadata['validation'] = report
    write_metadata(root / 'metadata.json', metadata)
    return root


def scientific_inputs(cfg, base):
    if cfg['force_field']['family'] not in ('gaff', 'gaff2') or cfg['force_field']['charge_model'] != 'resp':
        raise ValueError('Choose reviewed GAFF/GAFF2 + RESP; no AM1-BCC fallback')
    for key in ('qm_backend', 'method', 'basis', 'geometry_file', 'esp_conditions', 'fitting_conditions', 'source_file'):
        if not cfg['resp'].get(key): raise ValueError(f'Missing scientific input: resp.{key}')
    if cfg['resp']['source_kind'] != 'charges' or cfg['resp']['atom_order_confirmed'] is not True:
        raise ValueError('Prepare fitted RESP charges and explicitly confirm atom order')
    for key in ('geometry_file', 'source_file'):
        (base / cfg['resp'][key]).resolve(strict=True)
    nb = cfg['nonbonded']
    if nb['coulombtype'] != 'PME' or nb['vdwtype'] != 'Cut-off' or nb['dispersion_correction'] not in ('no', 'EnerPres', 'Ener'):
        raise ValueError('This adapter supports explicitly reviewed PME + Cut-off; other models require a separate adapter')
    for key in ('rcoulomb_nm', 'rvdw_nm', 'rlist_nm'):
        if not isinstance(nb[key], (int, float)) or not np.isfinite(nb[key]) or not 0 < nb[key] < .69:
            raise ValueError(f'{key} must be positive and below half the 1.38 nm box')
    if nb['rlist_nm'] < max(nb['rcoulomb_nm'], nb['rvdw_nm']) or not nb['adoption_reason']:
        raise ValueError('Require rlist >= cutoffs and a written adoption reason')
    if cfg['review']['force_field_parameters_confirmed'] is not True:
        raise ValueError('Review parameter/1-4 settings and confirm explicitly')
    threshold = cfg['review']['max_force_kJ_mol_nm']
    if not isinstance(threshold, (int, float)) or not np.isfinite(threshold) or threshold <= 0:
        raise ValueError('Provide an explicit reviewed force threshold')


def mdp(cfg, steps=600):
    nb = cfg['nonbonded']
    return f'''integrator = md-vv
 dt = 0.0005
 nsteps = {steps}
 constraints = none
 tcoupl = nose-hoover
 tc-grps = System
 tau-t = 0.3141592654
 ref-t = 300
 nh-chain-length = 3
 nsttcouple = 1
 pcoupl = no
 gen-vel = no
 continuation = no
 comm-mode = none
 nstcomm = 0
 pbc = xyz
 cutoff-scheme = Verlet
 nstlist = 1
 verlet-buffer-tolerance = -1
 rlist = {nb['rlist_nm']}
 coulombtype = {nb['coulombtype']}
 rcoulomb = {nb['rcoulomb_nm']}
 vdwtype = {nb['vdwtype']}
 rvdw = {nb['rvdw_nm']}
 DispCorr = {nb['dispersion_correction']}
 nstxout = {1 if steps == 0 else 0}
 nstvout = {1 if steps == 0 else 0}
 nstfout = {1 if steps == 0 else 0}
 nstxout-compressed = 5
 compressed-x-precision = 100000
 nstcalcenergy = 1
 nstenergy = 5
 nstlog = 5
'''


def prepare(cfg, config_path, build_dir, single_top, output, include_dirs=()):
    import parmed as pmd
    from mdworkbench.parameterization.resp import read_charges
    base = Path(config_path).resolve().parent
    scientific_inputs(cfg, base)
    order = cfg['structure']['ase_atom_order']
    if not isinstance(order, list) or sorted(order) != list(range(9)):
        raise ValueError('ase_atom_order must be a verified permutation of 0..8')
    source = Path(build_dir)
    data = np.load(source / 'initial.npz')
    report = json.loads((source / 'build_report.json').read_text())
    root = Path(output).resolve(); root.mkdir(parents=True, exist_ok=False)
    snapshot_topology(single_top, root / 'topology', include_dirs)
    with warnings.catch_warnings():
        warnings.simplefilter('error', pmd.exceptions.GromacsWarning)
        try:
            top = pmd.load_file(str(root / 'topology/system.top'))
        except pmd.exceptions.GromacsWarning as error:
            raise ValueError(f'Topology warning requires review: {error}') from error
    if len(top.atoms) != 9 or len(top.residues) != 1 or len(top.bonds) != 8:
        raise ValueError('Supply a single ethanol molecule topology (9 atoms, 8 bonds, one residue)')
    from ase.data import atomic_numbers
    if [a.atomic_number for a in top.atoms] != [atomic_numbers[report['ase_single_symbols'][i]] for i in order]:
        raise ValueError('Element order differs from declared ASE mapping')
    mapped_bonds = {tuple(sorted((order[b.atom1.idx], order[b.atom2.idx]))) for b in top.bonds}
    if mapped_bonds != {tuple(b) for b in report['ase_single_bonds']}:
        raise ValueError('Topology connectivity differs from ASE ethanol')
    charges = read_charges(base / cfg['resp']['source_file'], 9, 0)
    if not np.allclose([a.charge for a in top.atoms], np.asarray(charges)[order], atol=1e-5):
        raise ValueError('Topology charges do not match supplied RESP in ASE order')
    for terms in (top.bonds, top.angles, top.dihedrals, top.impropers):
        if any(t.type is None for t in terms): raise ValueError('Missing bonded parameter')
    text = (root/'topology/system.top').read_text()
    import re
    # Only the final [ molecules ] section can be changed; include graph is preserved.
    match = re.search(r'(?im)^\s*\[\s*molecules\s*\]\s*$', text)
    if not match: raise ValueError('Top-level topology must have [ molecules ]')
    rows = [line.split(';')[0].split() for line in text[match.end():].splitlines()]
    rows = [r for r in rows if r]
    if len(rows) != 1 or len(rows[0]) != 2 or rows[0][1] != '1':
        raise ValueError('Single topology must contain exactly one molecule count of 1')
    (root/'topology/system.top').write_text(text[:match.end()] + f'\n{rows[0][0]} 27\n')
    # Reject explicit constraint/virtual-site/rigid-water definitions anywhere in snapshot.
    for path in (root/'topology').iterdir():
        if re.search(r'\[\s*(constraints|settles|virtual_sites\w*)\s*\]', path.read_text(), re.I):
            raise ValueError('Topology contains constraints or virtual sites')
    permutation = np.concatenate([np.asarray(order) + 9*m for m in range(27)])
    xyz = data['positions_A'][permutation] / 10
    velocity = data['velocities_nm_ps'][permutation]
    if not np.allclose([a.mass for a in top.atoms], data['masses_amu'][order], atol=.02):
        raise ValueError('Topology masses differ materially from ASE masses')
    lines = ['Ethanol27: coordinates nm, velocities nm/ps', '243']
    for i, (pos, vel) in enumerate(zip(xyz, velocity)):
        atom = top.atoms[i % 9]
        lines.append(f'{i//9+1:5d}{"ETH":<5}{atom.name:>5}{i+1:5d}' + ''.join(f'{v:13.8f}' for v in pos) + ''.join(f'{v:13.9f}' for v in vel))
    lines.append(' 1.380000000 1.380000000 1.380000000')
    (root/'initial.gro').write_text('\n'.join(lines)+'\n')
    (root/'nvt.mdp').write_text(mdp(cfg))
    (root/'preflight.mdp').write_text(mdp(cfg, 0))
    shutil.copy2(config_path, root/'system.yaml')
    o = order.index(report['oxygen_ase_index']); h = order.index(report['hydroxyl_h_ase_index'])
    write_metadata(root/'atom_map.json', {'molecule_ids': np.repeat(np.arange(27),9).tolist(),
                   'oxygen': [9*m+o for m in range(27)], 'hydroxyl_h': [9*m+h for m in range(27)],
                   'ase_atom_order': order})
    ff = {'family': cfg['force_field']['family'],
          'defaults': {key: getattr(top.defaults, key) for key in ('nbfunc', 'comb_rule', 'gen_pairs', 'fudgeLJ', 'fudgeQQ')},
          'atoms': [{'name':a.name,'type':a.type,'charge':a.charge,'mass':a.mass,'sigma_A':a.sigma,'epsilon_kcal_mol':a.epsilon} for a in top.atoms]}
    for label, terms in [('bonds',top.bonds), ('angles',top.angles), ('dihedrals',top.dihedrals), ('impropers',top.impropers)]:
        ff[label] = [str(t) + ' type=' + str(t.type) for t in terms]
    ff['note'] = 'Full source parameters, exclusions, pairs and combination/1-4 rules: topology/*.top and *.itp; processed.top after grompp is authoritative.'
    write_metadata(root/'force_field_report.json', ff)
    inputs = [p for p in root.rglob('*') if p.is_file()] + [source/'initial.npz', base/cfg['resp']['source_file'], base/cfg['resp']['geometry_file']]
    meta = collect_metadata(inputs, 42, cfg)
    meta['build_report'] = report
    meta['snapshot_hashes'] = {str(p.relative_to(root)):hash_file(p) for p in root.rglob('*') if p.is_file()}
    meta['status'] = 'prepared'
    write_metadata(root/'metadata.json', meta)
    return root


def execute(root, args, input_text=None):
    from datetime import datetime, timezone
    from mdworkbench.gromacs.runner import detect_gromacs
    detected = detect_gromacs()
    if not detected.get('usable'): raise RuntimeError('GROMACS command unavailable; set GMX_BIN')
    command = [detected['executable'], *args]
    meta = json.loads((root/'metadata.json').read_text())
    entry = {'command':command,'cwd':str(root),'timestamp':datetime.now(timezone.utc).isoformat()}
    meta['gromacs'] = detected
    process = subprocess.run(command, cwd=root, input=input_text, text=True, capture_output=True)
    entry['returncode'] = process.returncode
    meta['commands'].append(entry)
    write_metadata(root/'metadata.json',meta)
    (root/(args[0] + '_' + str(len(meta['commands'])) + '.log')).write_text(process.stdout+process.stderr)
    if process.returncode: raise RuntimeError(f'Command failed; inspect {root}: {command}')


def verify_snapshot(root):
    meta = json.loads((root/'metadata.json').read_text())
    for name,digest in meta['snapshot_hashes'].items():
        if hash_file(root/name) != digest: raise ValueError(f'Prepared input changed: {name}; prepare a new run')
    return meta


def preflight(root):
    root = Path(root).resolve(); verify_snapshot(root)
    if (root/'preflight.tpr').exists(): raise ValueError('Preflight already attempted; prepare a new run')
    cfg = config(root/'system.yaml')
    execute(root, ['grompp','-f','preflight.mdp','-c','initial.gro','-p','topology/system.top','-o','preflight.tpr','-pp','processed.top','-po','preflight_processed.mdp'])
    execute(root, ['mdrun','-deffnm','preflight','-nt','1'])
    execute(root, ['energy','-f','preflight.edr','-o','preflight_energy.xvg'], 'Potential\n0\n')
    import MDAnalysis as mda
    u = mda.Universe(str(root/'initial.gro'),str(root/'preflight.trr'))
    ts = u.trajectory[0]
    force = np.linalg.norm(ts.forces,axis=1)*10 # MDAnalysis kJ/mol/A -> kJ/mol/nm
    energy = np.loadtxt(root/'preflight_energy.xvg',comments=['#','@'],ndmin=2)[:,1]
    mapping = json.loads((root/'atom_map.json').read_text())
    delta = ts.positions[:,None,:]-ts.positions[None,:,:]
    delta -= 13.8*np.rint(delta/13.8)
    distance = np.linalg.norm(delta,axis=-1)
    ids = np.asarray(mapping['molecule_ids'])
    closest = distance[ids[:,None]!=ids[None,:]].min()
    safe = bool(np.isfinite(force).all() and np.isfinite(energy).all() and force.max() <= cfg['review']['max_force_kJ_mol_nm'] and closest >= cfg['review']['catastrophic_distance_A'])
    hashes = {name:hash_file(root/name) for name in ('preflight.tpr','preflight.trr','preflight.edr','processed.top')}
    write_metadata(root/'preflight_report.json',{'safe':safe,'potential_kJ_mol':energy.tolist(),'max_force_kJ_mol_nm':float(force.max()),'closest_intermolecular_A':float(closest),'output_hashes':hashes})
    if not safe: raise ValueError('指定された初期条件ではGROMACS古典力場計算を安全に開始できない可能性がある。EMは実施しません。')


def run_nvt(root):
    root = Path(root).resolve(); verify_snapshot(root)
    report = json.loads((root/'preflight_report.json').read_text())
    if report['safe'] is not True: raise ValueError('Preflight must pass')
    for name,digest in report['output_hashes'].items():
        if hash_file(root/name)!=digest: raise ValueError('Preflight output changed')
    if (root/'nvt.tpr').exists(): raise ValueError('Run already attempted; prepare a new run')
    execute(root, ['grompp','-f','nvt.mdp','-c','initial.gro','-p','topology/system.top','-o','nvt.tpr','-pp','nvt_processed.top','-po','nvt_processed.mdp'])
    execute(root, ['mdrun','-deffnm','nvt','-nt','1'])
    execute(root, ['energy','-f','nvt.edr','-o','thermodynamics.xvg'], 'Temperature\nPotential\n0\n')
    meta=json.loads((root/'metadata.json').read_text()); meta['status']='completed'
    meta['output_hashes']={str(p.relative_to(root)):hash_file(p) for p in root.glob('nvt.*') if p.is_file()}
    write_metadata(root/'metadata.json',meta)


def analyze(root, mace=None, mace_steps=None):
    root=Path(root).resolve()
    mapping=json.loads((root/'atom_map.json').read_text())
    xyz,box,steps=read_gromacs(root/'nvt.tpr',root/'nvt.xtc')
    take=production_indices(steps)
    result=intermolecular_rdf(xyz[take],box[take],mapping['molecule_ids'],mapping['oxygen'],mapping['hydroxyl_h'])
    np.savez(root/'rdf_gromacs.npz',**result,steps=steps[take])
    if mace:
        xyz,box,steps=read_mace(mace, np.load(mace_steps) if mace_steps else None)
        take=production_indices(steps)
        # MACE reference must use molecule-major ASE atom order, unlike mapped GROMACS order.
        order=mapping['ase_atom_order']
        o=[9*m+order[mapping['oxygen'][0]] for m in range(27)]
        h=[9*m+order[mapping['hydroxyl_h'][0]] for m in range(27)]
        result=intermolecular_rdf(xyz[take],box[take],mapping['molecule_ids'],o,h)
        np.savez(root/'rdf_mace.npz',**result,steps=steps[take])
    return result


def read_thermodynamics(path):
    """Use XVG legends, since energy selection order need not equal column order."""
    import re
    text = Path(path).read_text()
    labels = {label: int(index) + 1 for index, label in
              re.findall(r'@\s+s(\d+)\s+legend\s+"([^"]+)"', text)}
    if not {'Temperature', 'Potential'} <= labels.keys():
        raise ValueError('XVG must contain labeled Temperature and Potential series')
    values = np.loadtxt(path, comments=['#', '@'], ndmin=2)
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite thermodynamics trace')
    return values[:, 0], values[:, labels['Temperature']], values[:, labels['Potential']]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['build','parameterize','prepare','preflight','run','analyze'])
    parser.add_argument('--config',default='systems/ethanol_27/system.yaml')
    parser.add_argument('--build-dir',default='runs/ethanol27_initial')
    parser.add_argument('--output')
    parser.add_argument('--single-top')
    parser.add_argument('--include-dir',action='append',default=[])
    parser.add_argument('--run-dir',default='runs/ethanol27_nvt')
    parser.add_argument('--mace'); parser.add_argument('--mace-steps')
    args=parser.parse_args()
    try:
        if args.action=='build': print(build(config(args.config),args.output or args.build_dir))
        elif args.action=='parameterize':
            from mdworkbench.config import load_config
            from mdworkbench.parameterization.amber import parameterize
            cfg=load_config(args.config)
            if not cfg['resp'].get('esp_conditions') or not cfg['resp'].get('fitting_conditions'):
                raise ValueError('Specify ESP and RESP fitting conditions before parameterization')
            print(parameterize(cfg,args.output or 'runs/ethanol27_parameters'))
        elif args.action=='prepare':
            if not args.single_top: raise ValueError('--single-top is required')
            print(prepare(config(args.config),args.config,args.build_dir,args.single_top,args.output or args.run_dir,args.include_dir))
        elif args.action=='preflight': preflight(args.run_dir)
        elif args.action=='run': run_nvt(args.run_dir)
        else: analyze(args.run_dir,args.mace,args.mace_steps)
    except (ValueError, RuntimeError, FileNotFoundError) as error:
        parser.exit(2,f'STOP: {error}\n')


if __name__=='__main__': main()
