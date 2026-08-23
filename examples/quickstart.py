#!/usr/bin/env python3
"""Self-contained quickstart — no external data needed.

Generates a 3D pose from a SMILES with RDKit, writes it to a PDB with *scrambled* atom names
(simulating a docking/co-folding output), then uses PoseMap to recover a specific atom by
chemical identity. Run:  python examples/quickstart.py
"""
import random
from rdkit import Chem
from rdkit.Chem import AllChem

from posemap import PoseMapper, parse_pdb


def make_scrambled_pdb(smiles, seed=0):
    """A PDB string whose ligand atom names are deliberately meaningless (like real tools)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    conf = mol.GetConformer()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    rng = random.Random(seed); rng.shuffle(heavy)
    lines, ctr = [], {}
    for serial, idx in enumerate(heavy, 1):
        a = mol.GetAtomWithIdx(idx); el = a.GetSymbol().upper()
        ctr[el] = ctr.get(el, 0) + rng.randint(1, 3)
        p = conf.GetAtomPosition(idx)
        lines.append(f"HETATM{serial:>5} {el+str(ctr[el]):<4} LIG C   1    "
                     f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2}")
    return "\n".join(lines) + "\nEND\n"


# 1-octen-3-one, an activated alkene (Michael acceptor). We want the beta-carbon.
smiles = "CCCCCC(=O)C=C"
pdb = make_scrambled_pdb(smiles)
print("The ligand's atom names in the 'pose' are meaningless:")
print("  " + "  ".join(ln[12:16].strip() for ln in pdb.splitlines() if ln.startswith("HETATM")))

struct = parse_pdb(pdb)
res = PoseMapper.from_smiles(smiles).map(struct)     # auto-detect the ligand

# The beta-carbon is the alkene C that is NOT attached to the carbonyl. Tag it :1 in SMARTS.
beta = res.atoms_by_smarts("[C:1]=[C][CX3]=[O]")     # returns a LIST (symmetry-aware)
print(f"\nmatched: {res.matched}   (isomorphisms: {len(res.mappings)})")
print(f"beta-carbon recovered: name={beta[0].name!r}  xyz={beta[0].xyz.round(3).tolist()}")

# Contrast: the carbonyl carbon
carb = res.atoms_by_smarts("[CX3:1]=[O]")
print(f"carbonyl carbon:       name={carb[0].name!r}  xyz={carb[0].xyz.round(3).tolist()}")
