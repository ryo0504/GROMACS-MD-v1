import pytest


def test_structure_generation_and_conversion(tmp_path):
    pytest.importorskip("rdkit")
    pytest.importorskip("ase")
    import numpy as np
    from mdworkbench.chemistry.build import build_molecule, write_sdf, read_molecule
    from mdworkbench.chemistry.convert import to_ase, update_coordinates
    a = build_molecule("CCO", 123)
    b = build_molecule("CCO", 123)
    np.testing.assert_allclose(a.GetConformer().GetPositions(), b.GetConformer().GetPositions())
    assert a.GetNumAtoms() == 9
    atoms = to_ase(a)
    assert atoms.get_chemical_formula() == "C2H6O"
    np.testing.assert_allclose(atoms.positions, a.GetConformer().GetPositions())
    atoms.positions += 1
    translated = update_coordinates(a, atoms)
    np.testing.assert_allclose(translated.GetConformer().GetPositions(), atoms.positions)
    file = tmp_path / "molecule.sdf"
    write_sdf(a, file)
    assert read_molecule(file).GetNumAtoms() == 9


def test_invalid_smiles():
    pytest.importorskip("rdkit")
    from mdworkbench.chemistry.build import build_molecule
    with pytest.raises(ValueError, match="SMILES"):
        build_molecule("not a molecule", 1)
