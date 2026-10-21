"""Explicit, single-conformer RHF/6-31G* → Amber RESP lesson for ethanol.

The surface is our deterministic Fibonacci-shell grid, not Gaussian's MK grid.
Only the isolated charge-fitting geometry is optimized. The MD initial state
is never edited. Completed stages are reused only after checking their hashes.
"""
from pathlib import Path
import json
import importlib.metadata
import numpy as np

from mdworkbench.chemistry.build import read_molecule
from mdworkbench.provenance.metadata import hash_file, write_metadata
from mdworkbench.parameterization.resp import read_charges
from mdworkbench.utils.subprocess import run


RADII_A = {"H": 1.20, "C": 1.70, "O": 1.52}
SCALES = (1.4, 1.6, 1.8, 2.0)


def ethanol_graph(sdf):
    """Require the explicit ethanol graph, preserving all nine atom identities."""
    mol = read_molecule(sdf)
    from rdkit import Chem
    if Chem.MolToSmiles(Chem.RemoveHs(mol)) != "CCO" or mol.GetNumAtoms() != 9:
        raise ValueError("This teaching recipe supports neutral, explicit-H ethanol only")
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    bonds = sorted(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds())
    groups = [[n.GetIdx() for n in a.GetNeighbors() if n.GetSymbol() == "H"]
              for a in mol.GetAtoms() if a.GetSymbol() == "C"]
    if sorted(map(len, groups)) != [2, 3]:
        raise ValueError("Expected CH3 and CH2 hydrogen groups")
    return mol, symbols, bonds, groups


def _cached(root, stage, signature):
    path = root / f"{stage}_report.json"
    if not path.exists():
        return None
    report = json.loads(path.read_text())
    if report["signature"] != signature:
        raise ValueError(f"{stage} inputs changed; choose a new output directory")
    for name, digest in report["output_hashes"].items():
        if hash_file(root / name) != digest:
            raise ValueError(f"{stage} output changed: {name}; choose a new directory")
    return report


def _save(root, stage, signature, files, **values):
    report = dict(signature=signature, output_hashes={name: hash_file(root / name) for name in files}, **values)
    write_metadata(root / f"{stage}_report.json", report)
    return report


def optimize(sdf, output):
    """Converged gas-phase RHF/6-31G* geometry, in unchanged SDF atom order."""
    import pyscf
    from pyscf import gto, scf
    from pyscf.geomopt.geometric_solver import kernel
    from rdkit import Chem

    sdf = Path(sdf).resolve()
    mol, symbols, bonds, groups = ethanol_graph(sdf)
    root = Path(output).resolve()
    signature = dict(sdf_sha256=hash_file(sdf), method="RHF", basis="6-31g*", charge=0,
                     spin=0, pyscf=pyscf.__version__, geometric=importlib.metadata.version("geometric"),
                     recipe_version=1)
    cached = _cached(root, "qm", signature)
    if cached:
        return cached
    root.mkdir(parents=True, exist_ok=False)
    (root / "input.sdf").write_bytes(sdf.read_bytes())
    xyz = mol.GetConformer().GetPositions()
    molecule = gto.M(atom=list(zip(symbols, xyz)), unit="Angstrom", basis="6-31g*",
                     charge=0, spin=0, symmetry=False, verbose=4,
                     output=str(root / "pyscf.log"), max_memory=2000)
    mf = scf.RHF(molecule)
    mf.conv_tol = 1e-10
    mf.max_cycle = 100
    mf.chkfile = str(root / "optimization.chk")
    converged, eq = kernel(mf, maxsteps=100, assert_convergence=True)
    if not converged:
        raise RuntimeError("Geometry did not converge; inspect pyscf.log")
    final = scf.RHF(eq)
    final.conv_tol = 1e-10
    final.max_cycle = 100
    final.chkfile = str(root / "hf.chk")
    energy = final.kernel()
    if not final.converged or not np.isfinite(energy):
        raise RuntimeError("Final SCF did not converge")
    coords = eq.atom_coords(unit="Angstrom")
    # Optimization must preserve the input bond graph (distance-based sanity check).
    from rdkit.Chem import rdDetermineBonds
    check = Chem.MolFromXYZBlock("9\noptimized ethanol\n" + "\n".join(
        f"{s} {x:.12f} {y:.12f} {z:.12f}" for s, (x, y, z) in zip(symbols, coords)))
    rdDetermineBonds.DetermineBonds(check, charge=0)
    found = sorted(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in check.GetBonds())
    if found != bonds:
        raise RuntimeError("Optimization changed the molecular connectivity")
    for i, pos in enumerate(coords):
        mol.GetConformer().SetAtomPosition(i, pos.tolist())
    with Chem.SDWriter(str(root / "optimized.sdf")) as writer:
        writer.write(mol)
    (root / "optimized.xyz").write_text("9\nRHF/6-31G*; original SDF atom order; Angstrom\n" + "\n".join(
        f"{s} {x:.12f} {y:.12f} {z:.12f}" for s, (x, y, z) in zip(symbols, coords)) + "\n")
    np.savez(root / "wavefunction.npz", coords_A=coords, density_matrix=final.make_rdm1())
    return _save(root, "qm", signature,
                 ["input.sdf", "optimized.sdf", "optimized.xyz", "wavefunction.npz", "hf.chk"],
                 geometry_converged=True, scf_converged=True, energy_hartree=float(energy),
                 symbols=symbols, bonds=bonds, hydrogen_groups=groups,
                 note="Single gas-phase conformer; no frequency/minimum certification; MD initial state unchanged")


def surface_points(coords_A, symbols, density=1.0):
    """Four exposed vdW shells; Fibonacci sphere, nominal points / Å²."""
    coords = np.asarray(coords_A, dtype=float)
    shells = []
    for scale in SCALES:
        radii = np.array([RADII_A[s] * scale for s in symbols])
        for i, radius in enumerate(radii):
            n = max(12, int(np.ceil(4 * np.pi * radius**2 * density)))
            k = np.arange(n)
            z = 1 - 2 * (k + .5) / n
            phi = k * np.pi * (3 - np.sqrt(5))
            r = np.sqrt(1 - z*z)
            points = coords[i] + radius * np.column_stack((r*np.cos(phi), r*np.sin(phi), z))
            distances = np.linalg.norm(points[:, None] - coords[None], axis=2)
            keep = np.all(distances >= radii[None] - 1e-10, axis=1)
            shells.append(points[keep])
    return np.concatenate(shells)


def potential(mol, dm, points_bohr):
    """V = Σ Z/r − ∫ rho(r')/|r-r'| dr', in atomic units.

    Uses the same Coulomb integrals as pyscf.tools.cubegen.mep.
    """
    from pyscf import gto, df
    points = np.asarray(points_bohr)
    distances = np.linalg.norm(points[:, None] - mol.atom_coords()[None], axis=2)
    nuclear = np.sum(mol.atom_charges()[None] / distances, axis=1)
    electronic = np.empty(len(points))
    for start in range(0, len(points), 256):
        stop = start + 256
        ints = df.incore.aux_e2(mol, gto.fakemol_for_charges(points[start:stop]))
        electronic[start:stop] = np.einsum("ijp,ij->p", ints, dm)
    return nuclear - electronic


def write_amber_esp(path, coords_bohr, grid_bohr, esp_au):
    """Amber espgen layout: 3I5 header, 17X3E16.7 nuclei, 4E16.7 ESP."""
    with Path(path).open("w") as handle:
        handle.write(f"{len(coords_bohr):5d}{len(grid_bohr):5d}{0:5d}\n")
        for x, y, z in coords_bohr:
            handle.write(f"{x:32.7E}{y:16.7E}{z:16.7E}\n")
        for v, (x, y, z) in zip(esp_au, grid_bohr):
            handle.write(f"{v:16.7E}{x:16.7E}{y:16.7E}{z:16.7E}\n")


def make_esp(output):
    from pyscf import gto
    from pyscf.data.nist import BOHR
    root = Path(output).resolve()
    qm = json.loads((root / "qm_report.json").read_text())
    _cached(root, "qm", qm["signature"])
    signature = dict(qm_report=hash_file(root / "qm_report.json"), radii_A=RADII_A,
                     scales=list(SCALES), density_per_A2=1.0, algorithm="Fibonacci exposed shells v1")
    cached = _cached(root, "esp", signature)
    if cached:
        return cached
    data = np.load(root / "wavefunction.npz")
    coords = data["coords_A"]
    molecule = gto.M(atom=list(zip(qm["symbols"], coords)), unit="Angstrom",
                     basis="6-31g*", charge=0, spin=0, symmetry=False, verbose=0)
    grid = surface_points(coords, qm["symbols"])
    esp = potential(molecule, data["density_matrix"], grid / BOHR)
    if len(grid) < 100 or not np.isfinite(esp).all():
        raise RuntimeError("Invalid ESP grid/potential")
    np.savez(root / "esp_grid.npz", grid_A=grid, esp_au=esp, coords_A=coords)
    write_amber_esp(root / "qm.esp", coords / BOHR, grid / BOHR, esp)
    return _save(root, "esp", signature, ["esp_grid.npz", "qm.esp"], point_count=len(grid),
                 units="coordinates: bohr in qm.esp, Angstrom in npz; potential: Hartree/e",
                 note="Custom exposed Fibonacci vdW shells; not identical to Gaussian MK sampling")


def resp_input(symbols, groups, stage):
    from rdkit.Chem import GetPeriodicTable
    flags = [0] * len(symbols)
    if stage == 2:
        flags = [-1 if s == "O" else 0 for s in symbols]
        carbon_h = {i for group in groups for i in group}
        for i, s in enumerate(symbols):
            if s == "H" and i not in carbon_h:
                flags[i] = -1  # keep hydroxyl H and O at stage-1 charges
        for group in groups:
            for i in sorted(group)[1:]:
                flags[i] = min(group) + 1  # RESP uses 1-based reference atom
    weight = .0005 if stage == 1 else .001
    text = (f"Ethanol single-conformer RESP stage {stage}\n\n &cntrl\n"
            f" nmol=1, ihfree=1, ioutopt=1, iqopt={1 if stage == 1 else 2},\n"
            f" qwt={weight},\n &end\n 1.0\nEthanol in original SDF atom order\n{0:5d}{9:5d}\n")
    for symbol, flag in zip(symbols, flags):
        text += f"{GetPeriodicTable().GetAtomicNumber(symbol):5d}{flag:5d}\n"
    return text + "\n\n"


def fit_resp(output):
    from pyscf.data.nist import BOHR
    root = Path(output).resolve()
    esp_report = json.loads((root / "esp_report.json").read_text())
    _cached(root, "esp", esp_report["signature"])
    _, symbols, _, groups = ethanol_graph(root / "input.sdf")
    signature = dict(esp_report=hash_file(root / "esp_report.json"), stages=[.0005, .001],
                     ihfree=1, hydrogen_groups=groups, stage2_fixed="O and hydroxyl H", recipe_version=1)
    cached = _cached(root, "resp", signature)
    if cached:
        return cached
    import shutil
    if not shutil.which("resp"):
        raise RuntimeError("AmberTools resp missing from PATH; use the gmx-md kernel")
    commands = []
    for stage in (1, 2):
        (root / f"resp{stage}.in").write_text(resp_input(symbols, groups, stage))
        command = ["resp", "-O", "-i", f"resp{stage}.in", "-o", f"resp{stage}.out",
                   "-e", "qm.esp", "-t", f"resp{stage}.qout", "-p", f"resp{stage}.pch"]
        if stage == 2:
            command += ["-q", "resp1.qout"]
        run(command, root, root / f"resp{stage}.log")
        commands.append(command)
        text = (root / f"resp{stage}.out").read_text()
        import re
        if not re.search(r"Convergence in\s+\d+ iterations", text) or "WARNING" in text.upper():
            raise RuntimeError(f"RESP stage {stage} requires review; inspect resp{stage}.out")
    charges = np.array(read_charges(root / "resp2.qout", 9, 0))
    first = np.array(read_charges(root / "resp1.qout", 9, 0))
    carbon_h = {i for group in groups for i in group}
    fixed = [i for i, s in enumerate(symbols) if s == "O" or (s == "H" and i not in carbon_h)]
    if not np.allclose(charges[fixed], first[fixed], atol=1e-6):
        raise RuntimeError("Stage-2 frozen charges changed")
    if any(np.ptp(charges[group]) > 1e-6 for group in groups):
        raise RuntimeError("Equivalent C-bound hydrogen charges differ")
    np.savetxt(root / "resp.charges", charges, fmt="%.10f")
    data = np.load(root / "esp_grid.npz")
    distance = np.linalg.norm(data["grid_A"][:, None] - data["coords_A"][None], axis=2) / BOHR
    predicted = (1 / distance) @ charges
    np.savez(root / "esp_fit.npz", qm_esp_au=data["esp_au"], fitted_esp_au=predicted)
    files = ["resp.charges", "esp_fit.npz"]
    files += [f"resp{s}.{suffix}" for s in (1, 2) for suffix in ("in", "out", "qout", "pch")]
    return _save(root, "resp", signature, files, charges=charges.tolist(), sum_charge=float(charges.sum()),
                 esp_rmse_au=float(np.sqrt(np.mean((predicted - data["esp_au"])**2))),
                 commands=commands, fixed_stage2_indices=fixed,
                 note="Two-stage educational single-conformer RESP; explicit custom ESP grid")


def parameter_config(template, output):
    """Write a new run-local config; never fill the generic system.yaml."""
    import yaml
    root = Path(output).resolve()
    qm = json.loads((root / "qm_report.json").read_text())
    fit = json.loads((root / "resp_report.json").read_text())
    _cached(root, "qm", qm["signature"])
    _cached(root, "resp", fit["signature"])
    cfg = yaml.safe_load(Path(template).read_text())
    cfg["force_field"]["family"] = "gaff2"
    cfg["resp"].update(qm_backend=f"PySCF {qm['signature']['pyscf']}", method="RHF",
        basis="6-31G*", geometry_file="optimized.xyz", source_kind="charges", source_file="resp.charges",
        esp_conditions=json.loads((root / "esp_report.json").read_text())["signature"],
        fitting_conditions=fit["signature"], atom_order_confirmed=True,
        fit_inputs=["resp1.in", "resp2.in"])
    cfg["structure"]["parameterization_file"] = "optimized.sdf"
    path = root / "parameterization.yaml"
    if path.exists() and yaml.safe_load(path.read_text()) != cfg:
        raise ValueError("Parameterization config changed; choose a new output directory")
    path.write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))
    return path


def parameterize(config_path, output):
    """Use the existing AmberTools pipeline; verify provenance on re-execution."""
    from mdworkbench.config import load_config
    from mdworkbench.parameterization.amber import parameterize as amber_parameterize
    root = Path(output).resolve()
    if not root.exists():
        return amber_parameterize(load_config(config_path), root)
    meta = json.loads((root / "metadata.json").read_text())
    if meta["scientific_config"] != load_config(config_path) or not meta.get("output_hashes"):
        raise ValueError("Incomplete/changed parameterization; choose a new output directory")
    for filename, digest in {**meta["input_hashes"], **meta["output_hashes"]}.items():
        if hash_file(filename) != digest:
            raise ValueError(f"Parameterization file changed: {filename}")
    return root


def review_topology(output, parameters):
    """Check coordinates as well as bonds/elements, including indistinguishable H."""
    import parmed as pmd
    from scipy.optimize import linear_sum_assignment
    from ase.data import atomic_numbers
    root, params = Path(output).resolve(), Path(parameters).resolve()
    candidates = list(params.rglob("*_GMX.top"))
    if len(candidates) != 1:
        raise ValueError("Expected exactly one ACPYPE single-molecule *_GMX.top")
    top_path = candidates[0]
    gro_path = top_path.with_suffix(".gro")
    top = pmd.load_file(str(top_path), xyz=str(gro_path))
    mol, symbols, bonds, groups = ethanol_graph(root / "optimized.sdf")
    coords = mol.GetConformer().GetPositions()
    if len(top.atoms) != 9 or len(top.residues) != 1:
        raise ValueError("Expected one nine-atom molecule")
    distances = np.linalg.norm(top.coordinates[:, None] - coords[None], axis=2)
    row, order = linear_sum_assignment(distances)
    if not np.array_equal(row, np.arange(9)) or distances[row, order].max() > .02:
        raise ValueError("Cannot establish coordinate-based topology→SDF atom mapping")
    if [a.atomic_number for a in top.atoms] != [atomic_numbers[symbols[i]] for i in order]:
        raise ValueError("Mapped elements differ")
    mapped = sorted(tuple(sorted((int(order[b.atom1.idx]), int(order[b.atom2.idx])))) for b in top.bonds)
    if mapped != bonds:
        raise ValueError("Mapped bond graph differs")
    charges = np.array(read_charges(root / "resp.charges", 9, 0))
    if not np.allclose([a.charge for a in top.atoms], charges[order], atol=1e-5):
        raise ValueError("Mapped RESP charges differ")
    charge_delta = np.array([a.charge for a in top.atoms]) - charges[order]
    if any(t.type is None for terms in (top.bonds, top.angles, top.dihedrals, top.impropers) for t in terms):
        raise ValueError("Missing bonded parameters")
    for atom in top.atoms:
        if not np.isfinite([atom.mass, atom.sigma, atom.epsilon, atom.charge]).all():
            raise ValueError("Nonfinite atomic parameter")
        if atom.mass <= 0 or atom.sigma < 0 or atom.epsilon < 0:
            raise ValueError("Invalid mass or LJ parameter")
    defaults = {k: getattr(top.defaults, k) for k in ("nbfunc", "comb_rule", "gen_pairs", "fudgeLJ", "fudgeQQ")}
    if (defaults["nbfunc"], defaults["comb_rule"], defaults["gen_pairs"]) != (1, 2, "yes"):
        raise ValueError("Unexpected GAFF2 combination/pair rules")
    if not np.isclose(defaults["fudgeLJ"], .5) or not np.isclose(defaults["fudgeQQ"], 1/1.2):
        raise ValueError("Unexpected GAFF2 1–4 scale factors")
    # Every graph-distance-three pair must appear explicitly; none may be added.
    graph = {i: set() for i in range(9)}
    for bond in top.bonds:
        i, j = bond.atom1.idx, bond.atom2.idx
        graph[i].add(j); graph[j].add(i)
    expected_pairs = set()
    for i in graph:
        visited, frontier = {i}, {i}
        for _ in range(3):
            frontier = set().union(*(graph[j] for j in frontier)) - visited
            visited |= frontier
        expected_pairs |= {(i, j) for j in frontier if i < j}
    actual_pairs = {tuple(sorted((t.atom1.idx, t.atom2.idx))) for t in top.adjusts}
    if actual_pairs != expected_pairs:
        raise ValueError("Explicit 1–4 pair list differs from molecular graph")
    import re
    match = re.search(r"\[\s*moleculetype\s*\]\s*\n(?:\s*;[^\n]*\n)*\s*\S+\s+(\d+)", top_path.read_text())
    if not match or int(match[1]) != 3:
        raise ValueError("Expected GAFF2 nrexcl=3 (1–2/1–3 excluded, 1–4 handled by pairs)")
    report = dict(checks_passed=True, single_top=str(top_path), ase_atom_order=order.tolist(), defaults=defaults,
        nrexcl=3, topology_charge_sum=float(sum(a.charge for a in top.atoms)),
        charge_rounding_delta_e=charge_delta.tolist(),
        coordinate_mapping_max_error_A=float(distances[row, order].max()),
        bonds=len(top.bonds), angles=len(top.angles), dihedrals=len(top.dihedrals), pairs_14=len(actual_pairs),
        atoms=[dict(topology_index=a.idx, ase_index=int(order[a.idx]), name=a.name, type=a.type,
                    charge=a.charge, mass=a.mass, sigma_A=a.sigma, epsilon_kcal_mol=a.epsilon) for a in top.atoms],
        source_hashes={str(p): hash_file(p) for p in (top_path, gro_path, root / "resp.charges")},
        note="Input consistency checks for this ethanol recipe; not validation of liquid properties")
    write_metadata(root / "topology_review.json", report)
    return report


def md_config(output, review):
    """Explicit small-box educational choices, after topology consistency checks."""
    import yaml
    root = Path(output).resolve()
    if review.get("checks_passed") is not True:
        raise ValueError("Topology checks must pass")
    for filename, digest in review["source_hashes"].items():
        if hash_file(filename) != digest:
            raise ValueError("Reviewed topology changed")
    cfg = yaml.safe_load((root / "parameterization.yaml").read_text())
    cfg["structure"]["ase_atom_order"] = review["ase_atom_order"]
    cfg["nonbonded"].update(coulombtype="PME", rcoulomb_nm=.60, vdwtype="Cut-off",
        rvdw_nm=.60, rlist_nm=.60, dispersion_correction="EnerPres",
        adoption_reason="操作学習用の1.38 nm固定箱で半箱0.69 nmを下回る0.60 nmを使用。"
        "PMEで長距離静電相互作用を扱い、LJに平均場の分散補正を適用。"
        "小cutoff・短時間・単一配座RESPであり、収束した液体物性の条件とは扱わない。")
    cfg["review"].update(force_field_parameters_confirmed=True, max_force_kJ_mol_nm=10000.)
    cfg["review"]["force_threshold_reason"] = (
        "操作学習用の発散予防screen。10000 kJ/mol/nmは1 uの原子にdt=0.5 fsで"
        "約0.0125 Åの半ステップ加速度変位を与える。安定性の保証ではなく、"
        "異常時は停止し、閾値引上げや初期系のEMは自動実施しない。")
    path = root / "system.yaml"
    if path.exists() and yaml.safe_load(path.read_text()) != cfg:
        raise ValueError("MD config changed; choose a new output directory")
    path.write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))
    return path
