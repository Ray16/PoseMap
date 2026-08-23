# PoseMap

**Recover chemically-meaningful atoms from a docking or co-folding pose — by identity, not by name.**

> Distribution name **PoseMap**; import as `posemap` (lowercase, per Python convention — like
> NumPy→`numpy`). CLI command: `posemap`.

Docking and co-folding tools (AutoDock/Gnina, **Boltz-1/2**, **AlphaFold3**, Chai, …) write ligand
atoms with *arbitrary* names and order, usually all under one generic residue name (`LIG`). So when
you need "the halo-carbon", "the enone β-carbon", "the anomeric carbon", or "the carbonyl oxygen"
out of a predicted structure, you can't just look it up by name — the name is meaningless.

`posemap` solves this by **graph isomorphism**: it maps a *template* molecule (a SMILES or RDKit mol,
where every atom's identity is known) onto the pose's ligand by matching **connectivity + element**.
You then ask for atoms with a **SMARTS** pattern on the template and get their **coordinates in the
pose** back — whatever the file happened to call them.

```python
from posemap import PoseMapper, load_structure

struct  = load_structure("pose.pdb")            # or .cif
mapper  = PoseMapper.from_smiles("CCCCCC(=O)C=C")
res     = mapper.map(struct, chain="C")         # or omit chain to auto-detect the ligand

beta = res.atoms_by_smarts("[C:1]=[C][CX3]=[O]") # Michael-acceptor β-carbon (:1)
print(beta[0].xyz)                               # -> array([-1.673, -0.341, -7.401])
```

Command line:

```bash
posemap pose.pdb "CCCCCC(=O)C=C" --motif michael_acceptor_beta --chain C
posemap pose.cif "O=C1C=CC(=O)N1" --smarts "[C:1]=[C][CX3]=[O]"   # symmetric -> 2 atoms
posemap pose.pdb "ClCC(=O)[O-]"   --dump                          # full template->pose table
posemap --list-motifs
```

## Why graph isomorphism (and not a distance heuristic)

The usual workaround is "the carbon nearest the halogen" / "the oxygen closest to the metal". Those
**drift onto the wrong atom** as poses vary — a whole class of silent bugs. `posemap` pins the atom by
**chemical identity**, so it is exact and pose-independent. Three properties fall out for free:

- **Name/order invariant** — works on any tool's scrambled ligand atoms.
- **Symmetry-aware** — a symmetric molecule (maleimide, benzene) returns the *set* of equivalent
  atoms; you pick the reactive orientation yourself (e.g. the copy closest to the catalytic residue).
- **A built-in pose gate** — if a pose's connectivity doesn't match the template (a clashed/distorted
  ligand), the mapping *fails cleanly* (`res.matched is False`) instead of returning a wrong atom.
  That failure is a useful signal: gate those poses out of your analysis.

## Design: two layers

1. **Generic mapper** (`PoseMapper`, `MapResult`) — chemistry-agnostic. Element-labelled VF2
   isomorphism (via `networkx`), covalent-radius bond perception, multi-ligand auto-detection by
   composition, symmetry as equivalent sets, connectivity gating. Reusable on *any* structure.
2. **Optional motif library** (`posemap.reactive`) — a thin, replaceable set of common reactive/
   functional-group SMARTS (`michael_acceptor_beta`, `halo_carbon`, `carbonyl_carbon`, …). Real
   projects keep their own mechanism-specific SMARTS; this layer just makes the tool useful out of
   the box and documents the pattern.

## API sketch

```python
mapper = PoseMapper.from_smiles(smiles)      # or .from_mol(rdkit_mol) for exact CCD bond orders
res    = mapper.map(structure, chain=None, resname=None, resseq=None, atoms=None)
res.matched          # bool — did the template map onto a ligand?
res.reason           # why not, if matched is False
res.mappings         # list of {template_idx -> pose_local_idx}; >1 == symmetry
res.atom(i)                       # -> [MappedAtom]  pose atoms for template atom i
res.atoms_by_smarts(smarts, 1)    # -> [MappedAtom]  pose atoms for the :1 SMARTS atom
res.candidates       # other composition-matching ligand copies present
mapper.map_all(structure)         # every matching ligand copy (e.g. multi-chain assemblies)

# MappedAtom: .xyz .name .element .local_idx .template_idx .serial
```

## Demonstrated at scale

`examples/demo_nac_systems.py` runs PoseMap over a whole enzyme-geometry project (13 systems, real
Boltz-2 co-folding poses): it auto-detects each substrate ligand by composition, maps it, and pulls
the mechanism-relevant reactive atom. Result — **97/97 substrates mapped, 100% of sampled poses**,
with symmetry sets and reactive atoms (halo-carbon, SNAr ipso-C, carbonyl C, enone β-C, thiopurine S)
recovered correctly. See `examples/EXAMPLE_OUTPUT.txt`.

## Install

```bash
pip install -e .          # needs numpy, networkx, rdkit
pip install -e '.[cif]'   # add gemmi for heavier mmCIF work (built-in reader covers common files)
pip install -e '.[test]' && pytest
```

## Scope & limits

- Handles **substrate/ligand** atoms. Protein-residue or standard-cofactor atoms (a catalytic Asp Oδ,
  FMN N5) keep standard names — select those by residue number / atom name directly; `posemap` is for
  the parts whose names you *can't* trust.
- Bond perception is distance + covalent-radii based; extreme distortions or unusual bonding may need a
  larger `bond_tol` or an RDKit-templated mol via `from_mol`.
- The built-in mmCIF reader parses the common `_atom_site` loop; for exotic CIF use `gemmi`.

MIT licensed. Contributions welcome — the motif library especially.
