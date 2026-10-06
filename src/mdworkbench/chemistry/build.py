from pathlib import Path


def build_molecule(smiles, seed):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError("RDKit conformer embedding failed")
    return mol


def read_molecule(path):
    from rdkit import Chem
    path = Path(path)
    if path.suffix.lower() == ".sdf":
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        if len(supplier) != 1:
            raise ValueError("v1 expects exactly one molecule in SDF")
        mol = supplier[0]
    elif path.suffix.lower() == ".mol2":
        mol = Chem.MolFromMol2File(str(path), removeHs=False)
    else:
        raise ValueError("Expected SDF or MOL2")
    if mol is None or not mol.GetNumConformers():
        raise ValueError("Structure must contain a valid explicit-H 3D conformer")
    if any(atom.GetNumImplicitHs() for atom in mol.GetAtoms()):
        raise ValueError("Input must contain explicit hydrogens; atom order must match RESP")
    return mol


def write_sdf(mol, path):
    from rdkit import Chem
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
