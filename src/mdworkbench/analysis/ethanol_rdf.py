"""Shared RDF: Angstrom coordinates, cubic box, directed intermolecular pairs."""
import numpy as np


def production_indices(steps):
    steps = np.asarray(steps)
    expected = np.arange(205, 601, 5)
    indices = np.flatnonzero(np.isin(steps, expected))
    if not np.array_equal(steps[indices], expected):
        raise ValueError('Require exactly steps 205,210,...,600 in increasing order (80 frames)')
    return indices


def intermolecular_rdf(coordinates, boxes, molecule_ids, oxygen, hydroxyl_h,
                       bins=45, r_max=6.762):
    xyz = np.asarray(coordinates, dtype=float)
    boxes = np.asarray(boxes, dtype=float)
    ids = np.asarray(molecule_ids)
    o, h = np.asarray(oxygen), np.asarray(hydroxyl_h)
    if xyz.ndim != 3 or xyz.shape[1:] != (len(ids), 3) or not len(xyz):
        raise ValueError('coordinates must be nonempty (frames, atoms, 3)')
    if boxes.shape != (len(xyz), 3) or not np.isfinite(xyz).all() or not np.isfinite(boxes).all():
        raise ValueError('Invalid/nonfinite coordinates or box')
    if not np.allclose(boxes, 13.8, atol=1e-4) or r_max > .5 * boxes.min():
        raise ValueError('This benchmark requires 13.8 A cubic boxes and r_max <= L/2')
    if len(o) != 27 or len(h) != 27 or len(set(ids[o])) != 27 or set(ids[o]) != set(ids[h]):
        raise ValueError('Require one oxygen and one hydroxyl H per each of 27 molecules')
    edges = np.linspace(0, r_max, bins + 1)
    shell = 4 * np.pi / 3 * np.diff(edges ** 3)
    result = {'r_A': (edges[1:] + edges[:-1]) / 2, 'edges_A': edges}
    for name, partner in [('OO', o), ('OH', h)]:
        mask = ids[o, None] != ids[partner][None, :]
        counts = np.zeros(bins)
        norm = np.zeros(bins)
        for frame, box in zip(xyz, boxes):
            delta = frame[o, None, :] - frame[partner][None, :, :]
            delta -= box * np.rint(delta / box)
            distances = np.linalg.norm(delta, axis=-1)[mask]
            counts += np.histogram(distances[distances < r_max], edges)[0]
            norm += 27 * (26 / np.prod(box)) * shell
        result[name] = counts / norm
        result[name + '_counts'] = counts
    return result


def read_gromacs(topology, trajectory):
    import MDAnalysis as mda
    u = mda.Universe(str(topology), str(trajectory))
    xyz, boxes, steps = [], [], []
    for ts in u.trajectory:
        step = int(round(ts.time / .0005))
        if abs(ts.time - step * .0005) > 1e-6:
            raise ValueError('Trajectory time is not on the 0.5 fs step grid')
        if not np.allclose(ts.dimensions[3:], 90):
            raise ValueError('Only orthogonal boxes supported')
        xyz.append(ts.positions.copy()); boxes.append(ts.dimensions[:3].copy()); steps.append(step)
    return np.asarray(xyz), np.asarray(boxes), np.asarray(steps)


def read_mace(path, saved_steps=None):
    from ase.io import read
    frames = read(str(path), ':')
    steps = saved_steps if saved_steps is not None else [a.info.get('step') for a in frames]
    if any(s is None for s in steps):
        raise ValueError('MACE frames need explicit step metadata or a saved_steps array; never infer stride')
    if any(not np.all(a.pbc) or not np.allclose(a.cell.array, np.eye(3) * 13.8) for a in frames):
        raise ValueError('MACE frames require cubic 13.8 A PBC')
    return np.asarray([a.positions for a in frames]), np.asarray([a.cell.lengths() for a in frames]), np.asarray(steps)
