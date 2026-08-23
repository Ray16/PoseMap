#!/usr/bin/env python3
"""Demonstration: PoseMap recovers substrate reactive atoms across a whole enzyme-geometry
project (the NAC_SN2 / preorg_to_activity systems), directly on real Boltz-2 co-folding poses.

For every modeled system it:
  1. auto-detects the substrate ligand in each pose by composition (no chain/name needed),
  2. maps the SMILES template onto it by graph isomorphism,
  3. reports the mapping success rate over a sample of poses (failures = distorted poses,
     correctly gated), and
  4. where a clean mechanism motif applies, pulls the reactive atom and prints its identity.

Run:  python examples/demo_nac_systems.py [SYSTEMS_DIR]
"""
import os, sys, csv, glob, random

from posemap import PoseMapper, load_structure
from posemap.reactive import motif_smarts

SYSTEMS = sys.argv[1] if len(sys.argv) > 1 else \
    "/nfs/lambda_stor_01/homes/rzhu/PMD/NAC_SN2/systems"

# system -> (reactive motif name or SMARTS, human label). None = map-only (motif is Layer-2 chemistry).
REACTIVE = {
    "oye":   ("[C:1]=[C][CX3]=[O]",   "hydride-acceptor beta-C"),
    "est2":  ("[CX3:1](=O)[OX2]",     "ester carbonyl C (Ser attack)"),
    "hld":   ("[CX4:1][Cl,Br,I]",     "halo-carbon (SN2 center)"),
    "ldex":  ("[CX4:1][Cl,Br,I]",     "halo-carbon (SN2 center)"),
    "gst":   ("[c:1][Cl,Br,I,F]",     "SNAr ipso-C (C-Cl)"),
    "tpmt":  ("[#16:1]",              "thiopurine S (methyl acceptor)"),
    "caffeine": (None, None),
    "galk":  (None, None), "smomt": (None, None), "pkinase": (None, None),
    "sult":  (None, None), "kemp": (None, None), "hheg": (None, None),
}
SAMPLE = 25   # poses sampled per system for the robustness rate


def smiles_lookup(csv_path):
    rows = list(csv.DictReader(open(csv_path)))
    if not rows:
        return {}, None
    col = "smiles" if "smiles" in rows[0] else "substrate_smiles"
    uniq = {r[col] for r in rows if r.get(col)}
    lut = {}
    for r in rows:
        smi = r.get(col)
        if not smi:
            continue
        for v in r.values():          # let any column value (member/substrate/enzyme) key the row
            if v:
                lut[v] = smi
    single = next(iter(uniq)) if len(uniq) == 1 else None
    return lut, single


def system_poses(sysdir):
    """Yield (outdir_name, pose_path) for every model_0 pose (one representative per substrate)."""
    for d in sorted(glob.glob(os.path.join(sysdir, "cofold", "out", "*"))):
        if not os.path.isdir(d):
            continue
        name = os.path.basename(d)
        ps = glob.glob(os.path.join(d, "boltz_results_*", "predictions", "*", "*model_0.pdb"))
        if ps:
            yield name, ps[0], os.path.dirname(ps[0])


def main():
    systems = [s for s in sorted(os.listdir(SYSTEMS))
               if os.path.isdir(os.path.join(SYSTEMS, s, "cofold", "out"))]
    print(f"PoseMap all-systems demonstration  (systems dir: {SYSTEMS})\n")
    hdr = f"{'system':10} {'subs':>5} {'mapped':>7} {'sym':>4} {'pose-rate':>10}  reactive-atom check"
    print(hdr); print("-" * len(hdr))
    rng = random.Random(0)
    grand_sub = grand_ok = 0

    for s in systems:
        sysdir = os.path.join(SYSTEMS, s)
        cand = glob.glob(os.path.join(sysdir, "data", "panel.csv")) + \
               glob.glob(os.path.join(sysdir, "data", "substrates.csv"))
        if not cand:
            continue
        lut, single = smiles_lookup(cand[0])
        reps = list(system_poses(sysdir))
        if not reps:
            continue

        n_sub = n_ok = n_sym = 0
        pose_pool = []          # (smiles, pose) for robustness sampling
        first_reactive = ""
        for name, pose, posedir in reps:
            smi = single or lut.get(name)
            if not smi:
                continue
            n_sub += 1
            try:
                mapper = PoseMapper.from_smiles(smi)
            except Exception as e:
                continue
            res = mapper.map(load_structure(pose))
            if res.matched:
                n_ok += 1
                if len(res.mappings) > 1:
                    n_sym += 1
                # reactive-atom check on the first mapped substrate of the system
                if not first_reactive and REACTIVE.get(s, (None,))[0]:
                    sm, label = REACTIVE[s]
                    hits = res.atoms_by_smarts(sm)
                    if hits:
                        h = hits[0]
                        extra = f" (+{len(hits)-1} sym)" if len(hits) > 1 else ""
                        first_reactive = f"{label}: {h.element} @[{h.xyz[0]:.1f},{h.xyz[1]:.1f},{h.xyz[2]:.1f}]{extra}"
                    else:
                        first_reactive = f"{label}: (SMARTS no match)"
            # collect a few poses for the rate
            all_poses = glob.glob(os.path.join(posedir, "*model*.pdb"))
            for p in all_poses:
                pose_pool.append((smi, p))

        # robustness: sample poses across the system, report mapping success
        rng.shuffle(pose_pool)
        sample = pose_pool[:SAMPLE]
        rate_ok = 0
        cache = {}
        for smi, p in sample:
            mp = cache.get(smi) or cache.setdefault(smi, PoseMapper.from_smiles(smi))
            if mp.map(load_structure(p)).matched:
                rate_ok += 1
        rate = f"{rate_ok}/{len(sample)}" if sample else "-"

        grand_sub += n_sub; grand_ok += n_ok
        print(f"{s:10} {n_sub:>5} {n_ok:>7} {n_sym:>4} {rate:>10}  {first_reactive}")

    print("-" * len(hdr))
    print(f"{'TOTAL':10} {grand_sub:>5} {grand_ok:>7}   substrates mapped "
          f"({100*grand_ok/max(grand_sub,1):.0f}% of representative poses)")


if __name__ == "__main__":
    main()
