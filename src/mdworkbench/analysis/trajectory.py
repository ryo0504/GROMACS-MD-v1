def load_trajectory(topology, trajectory):
    import MDAnalysis as mda
    return mda.Universe(str(topology), str(trajectory))


def rmsd(universe, selection="all"):
    """Make molecules whole and unwrap PBC beforehand for scientifically meaningful RMSD."""
    from MDAnalysis.analysis.rms import RMSD
    return RMSD(universe, select=selection, ref_frame=0).run().results.rmsd


def read_xvg(path):
    import numpy as np
    return np.loadtxt(path, comments=("#", "@"), ndmin=2)
