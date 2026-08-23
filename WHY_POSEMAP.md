# Why PoseMap? (vs. "just ask an agent" or "click it in PyMOL")

In the agentic era the fair question is: *why does this package need to exist — can't an LLM agent
just read the structure, or a human pick the atom in PyMOL?* Here is the honest, measured answer.

## The task

Get one specific **chemical** atom — e.g. the hydride-acceptor **β-carbon** of an ene-reductase
substrate — out of a docking / co-folding pose. The catch: tools like Boltz-2 and AlphaFold3 assign
ligand atom **names arbitrarily** (`C16`, `C17`, …) and drop all bond-order information. So "the
β-carbon" is defined by *chemistry*, but the file only gives you *coordinates and meaningless labels*.

## The evidence (measured on 400 real Boltz-2 OYE poses — `examples/why_posemap.py`)

| # | Question | Result |
|---|----------|--------|
| 1 | Can you select the β-carbon by atom **name**? | **No.** Across 8 substrates the β-carbon appears under **8 different names** (C8, C9, C10, C12, C14, C15, C17, C18). No convention exists. |
| 2 | Is the plausible **heuristic an agent would write** ("β-C = substrate carbon nearest the flavin N5") reliable? | **No — silently wrong on 34% of poses** (138/400), up to **88–96%** for some substrates. No error is ever raised. |
| 3 | Does the naive pick handle **symmetry**? | **No.** For maleimide (two equivalent β-carbons) it flips between them pose-to-pose and sometimes lands on a **carbonyl** carbon. PoseMap returns both equivalents on 50/50 poses. |
| 4 | Does it **scale** and reproduce? | PoseMap: 400 poses in ~29 s, deterministic, bit-for-bit reproducible. |

The 34% figure is the crux: the wrong answers are **plausible and silent** — exactly the bugs that
survive review and corrupt a downstream correlation. (This is the same failure class as the "nearest
carboxylate drifted off the catalytic Asp" selector artifact that once produced a spurious ρ=−1.00 in
a sister project.)

## "Why not just ask an LLM agent?"

An agent can attempt this three ways, and all three are worse than a library call:

1. **Read the file and reason.** Atom names carry no information, so the agent must reconstruct
   connectivity from raw coordinates — a 3D geometric-reasoning task LLMs are unreliable at, and one
   that yields **no signal when it's wrong**. (We ran this experiment directly; see below.)
2. **Write an ad-hoc heuristic** (nearest-atom, name match). That is precisely the "nearest N5 carbon"
   rule above: **silently wrong 34% of the time here**, and it degrades unpredictably per substrate.
3. **Do it correctly** — which means implementing element-labelled graph isomorphism + a SMARTS query
   + a connectivity gate… i.e. **re-deriving PoseMap on every call**, non-deterministically and at
   token cost. The correct agent move is therefore to *call* PoseMap, not to reinvent it each time.

Two more properties an agent can't give you for free:
- **Determinism / reproducibility.** Science needs the same atom every run; an LLM call does not.
- **A loud failure mode.** PoseMap returns `matched is False` when a pose's connectivity doesn't match
  the template (distorted/wrong pose) — a built-in correctness gate. An agent fails *silently* with a
  confident-looking wrong coordinate.

### Direct agent experiment (honest result)

We gave a capable LLM agent the raw coordinates of a simple substrate (1-octen-3-one) and a symmetric
one (maleimide) and asked for the β-carbon, **forbidding cheminformatics libraries** — the "just read
the structure and reason" path a user would try.

**It got both right**, including recognizing that maleimide has two symmetry-equivalent β-carbons
(C8≡C9). That is the honest outcome — and it is exactly the point. It succeeded *only by manually
re-deriving PoseMap's algorithm*: it computed interatomic distances, inferred bond orders (C=O ≈ 1.22 Å,
C=C ≈ 1.34–1.41 Å, C–C ≈ 1.45–1.54 Å), traced the connectivity, and identified each atom by its
chemical role — a ~14k-token, ~20-second reasoning episode **per molecule**, non-deterministic, with
**no correctness gate**.

So the realistic picture at scale is a three-way choice:

| Approach | Correct? | Cost / pose | Deterministic? | Fails loudly on a bad pose? |
|---|---|---|---|---|
| Select by atom **name** | impossible (Exp 1) | — | — | — |
| **Nearest-atom** heuristic | **wrong 34%** (Exp 2) | ~0 | yes | no (silent) |
| Agent **reasons from coordinates** | right here, but re-derives the algorithm every call | ~14k tokens, ~20 s | **no** | no (confident even when wrong) |
| **PoseMap** | right (identity + gate) | ~**70 ms** | **yes** | **yes** (`matched is False`) |

The agent's *correct* path **is** PoseMap's algorithm, re-implemented from scratch on every call.
PoseMap makes that one deterministic, gated, ~300× faster call — and on a distorted pose it returns
`matched is False` instead of a confident wrong coordinate.

## "Why not pick the ligand/atom manually in PyMOL?"

PyMOL is excellent for *looking* at one structure. It is the wrong tool for *extracting a defined atom
programmatically, at scale*:

- **Scale.** 400 poses × dozens of substrates × many enzymes = tens of thousands of picks. Manual
  clicking is minutes each; the analysis above would take weeks by hand and seconds with PoseMap.
- **Reproducibility & provenance.** A pipeline (or a paper's methods) can't contain "a human clicked
  the atom." PoseMap's selection is a versioned SMARTS string — auditable and re-runnable.
- **Symmetry consistency.** A human clicking maleimide picks *one* of the two equivalent β-carbons,
  differently on different poses. PoseMap returns the equivalence set explicitly so you resolve it by
  a stated rule, not by luck.
- **Novel scaffolds.** Eyeballing "which carbon is β" on an unfamiliar molecule is error-prone; a
  SMARTS pattern is an unambiguous definition that travels to any new substrate.
- **It defeats automation.** The entire point of an autonomous enzyme-analysis agent is that no human
  is in the per-pose loop.

## Honest scope (what PoseMap is *not* claiming)

PoseMap does nothing that RDKit's `AssignBondOrdersFromTemplate` / substructure matching, or `spyrmsd`,
could not in principle do. Its value is **packaging**: a robust, connectivity-gated, symmetry-aware,
tool- and format-agnostic, one-call implementation — plus a reactive-motif library and CLI — so that
pipelines and agents get the atom **right, deterministically, without reinventing it or subtly getting
it wrong**. In the agentic era, that reliability *is* the product: it's the difference between an agent
that re-derives a fiddly algorithm (badly, 34% of the time) on every call and one that makes a single
correct, free, reproducible call.
