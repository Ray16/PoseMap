"""Optional convenience layer: a small library of common reactive-/functional-atom SMARTS.

This is deliberately thin and separate from the generic mapper. The mapper (``PoseMapper`` +
``MapResult.atoms_by_smarts``) is chemistry-agnostic; real projects should keep their own
mechanism-specific SMARTS. These are provided so ``posemap`` is useful out of the box and to
document the intended pattern: each entry maps a chemically-meaningful atom via map number ``1``.
"""
from __future__ import annotations

# name -> (SMARTS with a :1-tagged target atom, human description)
MOTIFS = {
    "michael_acceptor_beta": ("[C:1]=[C][CX3]=[O]", "beta-carbon of an enone/enal (hydride/Nu acceptor)"),
    "michael_acceptor_alpha": ("[C]=[C:1][CX3]=[O]", "alpha-carbon of an enone/enal"),
    "carbonyl_carbon": ("[CX3:1]=[OX1]", "carbonyl carbon (Buergi-Dunitz electrophile)"),
    "carbonyl_oxygen": ("[CX3]=[OX1:1]", "carbonyl oxygen"),
    "halo_carbon": ("[CX4:1][F,Cl,Br,I]", "sp3 carbon bearing a halogen (SN2 electrophile)"),
    "halide": ("[CX4][F,Cl,Br,I:1]", "halide leaving group"),
    "hydroxyl_oxygen": ("[OX2H:1]", "hydroxyl oxygen (nucleophile)"),
    "carboxylate_oxygen": ("[OX1:1][CX3]=[OX1,OX2]", "carboxylate/carboxyl oxygen"),
    "primary_amine_n": ("[NX3;H2:1]", "primary amine nitrogen"),
    "thiol_sulfur": ("[SX2H:1]", "thiol sulfur (nucleophile)"),
    "anomeric_carbon": ("[C:1]([OX2])[OX2]", "anomeric carbon (two O substituents)"),
    "aromatic_halide_ipso": ("[c:1][F,Cl,Br,I]", "ipso aromatic carbon bearing a halide (SNAr)"),
}


def motif_smarts(name: str) -> str:
    if name not in MOTIFS:
        raise KeyError(f"unknown motif {name!r}; known: {', '.join(sorted(MOTIFS))}")
    return MOTIFS[name][0]


def find_motif(result, name: str):
    """Convenience: ``result.atoms_by_smarts(motif_smarts(name))``."""
    return result.atoms_by_smarts(motif_smarts(name), mapnum=1)
