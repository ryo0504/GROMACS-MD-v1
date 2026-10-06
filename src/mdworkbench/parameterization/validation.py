from pathlib import Path
import math


def validate_frcmod(path):
    text = Path(path).read_text()
    if not text.strip() or "ATTN" in text.upper() or "NEED REVISION" in text.upper():
        raise ValueError("parmchk2 missing/estimated parameter requires review; inspect frcmod")


def topology_summary(topology, coordinates=None, net_charge=None, atom_count=None):
    import parmed as pmd
    structure = pmd.load_file(str(topology), xyz=str(coordinates)) if coordinates else pmd.load_file(str(topology))
    if not structure.atoms or any(not atom.type for atom in structure.atoms):
        raise ValueError("Empty topology or missing atom type")
    charges = [float(a.charge) for a in structure.atoms]
    if not all(math.isfinite(q) for q in charges):
        raise ValueError("Nonfinite charge")
    if atom_count is not None and len(structure.atoms) != atom_count:
        raise ValueError("Topology/structure atom count mismatch")
    if net_charge is not None and abs(sum(charges) - net_charge) > 1e-4:
        raise ValueError("Topology net charge mismatch")
    if any(b.type is None for b in structure.bonds) or any(a.type is None for a in structure.angles):
        raise ValueError("Missing bonded parameters")
    return {"atom_count": len(structure.atoms), "residue_count": len(structure.residues),
            "total_charge": sum(charges), "atoms": [{"name": a.name, "type": a.type,
            "charge": a.charge} for a in structure.atoms]}
