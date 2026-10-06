def to_ase(mol):
    from ase import Atoms
    return Atoms(symbols=[a.GetSymbol() for a in mol.GetAtoms()],
                 positions=mol.GetConformer().GetPositions())


def update_coordinates(mol, atoms):
    from rdkit.Geometry import Point3D
    if [a.GetSymbol() for a in mol.GetAtoms()] != atoms.get_chemical_symbols():
        raise ValueError("Atom identity/order changed")
    copy = type(mol)(mol)
    for index, xyz in enumerate(atoms.positions):
        copy.GetConformer().SetAtomPosition(index, Point3D(*map(float, xyz)))
    return copy
