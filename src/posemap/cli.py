"""Command-line interface: map a SMILES onto a structure and print the atoms you ask for.

    posemap POSE.pdb "CCCCCC(=O)C=C" --smarts "[C:1]=[C][CX3]=[O]"
    posemap POSE.cif "O=C1C=CC(=O)N1" --motif michael_acceptor_beta --chain C
    posemap POSE.pdb "..." --dump          # print the full template->pose atom table
"""
from __future__ import annotations

import argparse
import sys

from .io import load_structure
from .mapper import PoseMapper
from .reactive import motif_smarts, MOTIFS


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="posemap", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("structure", help="PDB/mmCIF file")
    p.add_argument("smiles", help="template SMILES")
    p.add_argument("--smarts", help="SMARTS with a :1-tagged target atom")
    p.add_argument("--motif", choices=sorted(MOTIFS), help="named motif instead of --smarts")
    p.add_argument("--mapnum", type=int, default=1, help="SMARTS atom-map number to return")
    p.add_argument("--chain", help="restrict to this chain")
    p.add_argument("--resname", help="restrict to this residue name")
    p.add_argument("--resseq", help="restrict to this residue number")
    p.add_argument("--bond-tol", type=float, default=0.45, help="covalent bond tolerance (Å)")
    p.add_argument("--dump", action="store_true", help="print the full template->pose table")
    p.add_argument("--list-motifs", action="store_true", help="list built-in motifs and exit")
    args = p.parse_args(argv)

    if args.list_motifs:
        for name, (sm, desc) in sorted(MOTIFS.items()):
            print(f"{name:26} {sm:26} {desc}")
        return 0

    struct = load_structure(args.structure)
    mapper = PoseMapper.from_smiles(args.smiles, bond_tol=args.bond_tol)
    res = mapper.map(struct, chain=args.chain, resname=args.resname, resseq=args.resseq)
    if not res.matched:
        print(f"NO MATCH: {res.reason}", file=sys.stderr)
        return 2

    print(f"# matched ligand {res.ligand_key}  ({len(res.mappings)} isomorphism(s); "
          f"{'symmetric' if len(res.mappings) > 1 else 'unique'})")
    if res.candidates:
        print(f"# note: other composition-matching ligands present: {res.candidates}")

    if args.dump:
        print(f"{'tmpl_idx':>8} {'element':>7} {'pose_name':>9} {'x':>8} {'y':>8} {'z':>8}")
        for ti in range(res.template.GetNumAtoms()):
            for ma in res.atom(ti):
                print(f"{ti:>8} {ma.element:>7} {ma.name:>9} "
                      f"{ma.xyz[0]:>8.3f} {ma.xyz[1]:>8.3f} {ma.xyz[2]:>8.3f}")
        return 0

    smarts = args.smarts or (motif_smarts(args.motif) if args.motif else None)
    if not smarts:
        print("provide --smarts, --motif, or --dump", file=sys.stderr)
        return 1
    hits = res.atoms_by_smarts(smarts, mapnum=args.mapnum)
    if not hits:
        print(f"# SMARTS {smarts!r} matched no atoms", file=sys.stderr)
        return 3
    print(f"{'pose_name':>9} {'element':>7} {'x':>8} {'y':>8} {'z':>8}")
    for ma in hits:
        print(f"{ma.name:>9} {ma.element:>7} {ma.xyz[0]:>8.3f} {ma.xyz[1]:>8.3f} {ma.xyz[2]:>8.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
