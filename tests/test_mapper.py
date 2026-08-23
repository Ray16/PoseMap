"""Self-contained tests: generate a 3D pose from a SMILES with RDKit, write it out with
DELIBERATELY SCRAMBLED atom names/order (simulating what a docking/co-folding tool emits),
then assert posemap recovers the correct chemical atom by identity."""
import random

import numpy as np
import pytest

from rdkit import Chem
from rdkit.Chem import AllChem

from posemap import PoseMapper, parse_pdb, parse_cif
from posemap.reactive import motif_smarts


def _scrambled_structure(smiles, seed=0, fmt="pdb", resname="LIG", chain="C"):
    """Return (Structure, truth) where truth[serial] = original RDKit heavy-atom idx.
    Heavy atoms are written in shuffled order with arbitrary names Cxx/Nxx/Oxx."""
    rng = random.Random(seed)
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    conf = mol.GetConformer()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    order = heavy[:]
    rng.shuffle(order)
    truth = {}
    lines = []
    counters = {}
    for serial, idx in enumerate(order, start=1):
        a = mol.GetAtomWithIdx(idx)
        el = a.GetSymbol().upper()
        counters[el] = counters.get(el, 0) + rng.randint(1, 3)  # arbitrary, non-sequential
        name = f"{el}{counters[el]}"
        p = conf.GetAtomPosition(idx)
        truth[serial] = idx
        if fmt == "pdb":
            # strict PDB columns: name 13-16, altLoc 17, resName 18-20, chain 22, resSeq 23-26
            lines.append(
                f"HETATM{serial:>5} {name:<4} {resname:>3} {chain}{1:>4}    "
                f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2}")
        else:  # mmcif atom_site row
            lines.append(f"HETATM {serial} {el} {name} . {resname} {chain} 1 {1} "
                         f"{p.x:.3f} {p.y:.3f} {p.z:.3f} 1.00 0.00 1")
    if fmt == "pdb":
        text = "\n".join(lines) + "\nEND\n"
        return parse_pdb(text), truth, mol
    header = ("data_test\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
              "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"
              "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.label_seq_id\n"
              "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
              "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_PDB_model_num\n")
    return parse_cif(header + "\n".join(lines) + "\n#\n"), truth, mol


def _truth_indices(mol, smarts, mapnum=1):
    q = Chem.MolFromSmarts(smarts)
    qpos = [a.GetIdx() for a in q.GetAtoms() if a.GetAtomMapNum() == mapnum][0]
    return {m[qpos] for m in mol.GetSubstructMatches(q)}


def test_unique_beta_carbon_pdb():
    smiles = "CCCCCC(=O)C=C"        # 1-octen-3-one; beta is the terminal =CH2
    struct, truth, mol = _scrambled_structure(smiles, seed=1)
    res = PoseMapper.from_smiles(smiles).map(struct, chain="C")
    assert res.matched and len(res.mappings) == 1
    hits = res.atoms_by_smarts(motif_smarts("michael_acceptor_beta"))
    assert len(hits) == 1
    # the recovered pose atom must correspond to a true beta atom of the template
    truth_beta = _truth_indices(mol, motif_smarts("michael_acceptor_beta"))
    assert truth[hits[0].serial] in truth_beta


def test_symmetry_returns_equivalent_set():
    smiles = "O=C1C=CC(=O)N1"        # maleimide: the two alkene carbons are equivalent
    struct, truth, mol = _scrambled_structure(smiles, seed=3)
    res = PoseMapper.from_smiles(smiles).map(struct, chain="C")
    assert res.matched
    hits = res.atoms_by_smarts(motif_smarts("michael_acceptor_beta"))
    assert len(hits) == 2                       # both symmetry-equivalent beta carbons
    truth_beta = _truth_indices(mol, motif_smarts("michael_acceptor_beta"))
    assert {truth[h.serial] for h in hits} == truth_beta


def test_halo_carbon_and_halide():
    smiles = "ClCC(=O)[O-]"          # chloroacetate (L-DEX substrate)
    struct, truth, mol = _scrambled_structure(smiles, seed=5)
    res = PoseMapper.from_smiles(smiles).map(struct, chain="C")
    assert res.matched
    c = res.atoms_by_smarts(motif_smarts("halo_carbon"))
    x = res.atoms_by_smarts(motif_smarts("halide"))
    assert len(c) == 1 and c[0].element == "C"
    assert len(x) == 1 and x[0].element == "CL"


def test_gate_rejects_broken_connectivity():
    smiles = "CCCCCC(=O)C=C"
    struct, truth, mol = _scrambled_structure(smiles, seed=1)
    # tear a bond: push the first heavy atom 5 Å away so its connectivity changes
    struct.atoms[0].xyz = struct.atoms[0].xyz + np.array([5.0, 5.0, 5.0])
    res = PoseMapper.from_smiles(smiles).map(struct, chain="C")
    assert not res.matched
    assert "isomorphism" in res.reason or "count" in res.reason


def test_auto_detection_without_selection():
    smiles = "O=C1CCC=C1"            # 2-cyclopenten-1-one
    struct, truth, mol = _scrambled_structure(smiles, seed=7)
    res = PoseMapper.from_smiles(smiles).map(struct)   # no chain/resname given
    assert res.matched
    assert res.atoms_by_smarts(motif_smarts("carbonyl_carbon"))


def test_cif_roundtrip():
    smiles = "O=C1CCC=C1"
    struct, truth, mol = _scrambled_structure(smiles, seed=7, fmt="cif")
    assert len(struct) > 0
    res = PoseMapper.from_smiles(smiles).map(struct, chain="C")
    assert res.matched
