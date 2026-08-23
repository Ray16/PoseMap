"""posemap — recover chemically-meaningful atoms from a structure by graph isomorphism.

Docking and co-folding tools (AutoDock/Gnina, Boltz, AlphaFold3, ...) emit ligand atoms with
arbitrary names and order, all under a generic residue name. ``posemap`` maps a template
molecule (a SMILES or RDKit mol, where you *know* which atom is which) onto that structure by
element-labelled graph isomorphism, so you can ask for "the halo-carbon", "the enone
beta-carbon", "the anomeric carbon" and get its coordinate back — regardless of what the file
called it. Poses whose connectivity does not match the template fail cleanly (a gate for
distorted poses), and symmetric molecules return the equivalent atom set.

Quick start
-----------
    from posemap import PoseMapper, load_structure

    struct = load_structure("pose.pdb")
    mapper = PoseMapper.from_smiles("CCCCCC(=O)C=C")   # 1-octen-3-one
    res = mapper.map(struct, chain="C")                # or auto-detect the ligand
    beta = res.atoms_by_smarts("[C:1]=[C][CX3]=[O]")   # -> [MappedAtom(...)]
    print(beta[0].xyz)
"""
from .io import Structure, Atom, load_structure, parse_pdb, parse_cif
from .mapper import PoseMapper, MapResult, MappedAtom
from .graph import build_graph, covalent_radius
from . import reactive

__version__ = "0.1.0"

__all__ = [
    "Structure", "Atom", "load_structure", "parse_pdb", "parse_cif",
    "PoseMapper", "MapResult", "MappedAtom",
    "build_graph", "covalent_radius", "reactive",
    "map_smiles_to_structure",
]


def map_smiles_to_structure(smiles: str, path_or_structure, **map_kwargs) -> MapResult:
    """One-shot convenience: ``PoseMapper.from_smiles(smiles).map(structure, **kwargs)``."""
    struct = (path_or_structure if isinstance(path_or_structure, Structure)
              else load_structure(path_or_structure))
    return PoseMapper.from_smiles(smiles).map(struct, **map_kwargs)
