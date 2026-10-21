"""Independent electrostatics and file-boundary checks for the local QM lesson."""
from pathlib import Path
import numpy as np
import pytest

from mdworkbench.parameterization.pyscf_ethanol import potential, surface_points, write_amber_esp, _cached, _save


def test_potential_matches_analytic_gaussian_electron():
    """One normalized s Gaussian: electron potential = erf(sqrt(2α)r)/r."""
    pytest.importorskip("pyscf")
    from pyscf import gto
    from scipy.special import erfc
    atom = gto.M(atom="H 0 0 0", basis={"H": [[0, [1.0, 1.0]]]}, spin=1, verbose=0)
    radii = np.array([.4, 1.0, 2.0, 5.0])  # bohr, not Angstrom
    points = np.column_stack((radii, np.zeros((len(radii), 2))))
    actual = potential(atom, np.ones((1, 1)), points)
    expected = erfc(np.sqrt(2) * radii) / radii
    assert np.allclose(actual, expected, atol=1e-12)


def test_surface_is_exposed_and_amber_esp_preserves_units(tmp_path):
    coords = np.array([[0., 0., 0.], [2., 0., 0.]])
    points = surface_points(coords, ["C", "O"])
    # No surface point may be inside either atom's smallest shell.
    distances = np.linalg.norm(points[:, None] - coords[None], axis=2)
    assert np.all(distances >= np.array([1.7, 1.52]) * 1.4 - 1e-10)
    file = tmp_path / "qm.esp"
    write_amber_esp(file, [[1., -2., 3.]], [[4., 5., -6.]], [-.0123])
    lines = file.read_text().splitlines()
    assert [int(lines[0][i:i+5]) for i in (0, 5, 10)] == [1, 1, 0]
    # Read the same fixed columns as Amber, independent of the writer logic.
    assert lines[1][:16].isspace()
    assert [float(lines[1][i:i+16]) for i in (16, 32, 48)] == [1., -2., 3.]
    assert [float(lines[2][i:i+16]) for i in (0, 16, 32, 48)] == [-.0123, 4., 5., -6.]


def test_cached_scientific_result_rejects_modified_data(tmp_path):
    (tmp_path / "charges").write_text("0.1 -0.1\n")
    _save(tmp_path, "test", {"method": "test-only"}, ["charges"])
    assert _cached(tmp_path, "test", {"method": "test-only"})
    with pytest.raises(ValueError, match="inputs changed"):
        _cached(tmp_path, "test", {"method": "different"})
    (tmp_path / "charges").write_text("0 0\n")
    with pytest.raises(ValueError, match="output changed"):
        _cached(tmp_path, "test", {"method": "test-only"})
