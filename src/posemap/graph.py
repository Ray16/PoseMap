"""Covalent-bond perception and molecular-graph construction.

Bonds are inferred from interatomic distance against summed covalent radii with a tolerance
(the standard ``d < r_i + r_j + tol`` rule). This is more robust across chemistries than a
single fixed cutoff: it bonds a long C-I (~2.14 Å) without spuriously bonding a 1.9 Å
non-bonded contact between two carbons.
"""
from __future__ import annotations

import numpy as np
import networkx as nx

# Covalent radii (Å), Cordero et al. 2008 (single-bond) — organics + common metals.
COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "SI": 1.11,
    "P": 1.07, "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39, "SE": 1.20, "AS": 1.19,
    "TE": 1.38, "NA": 1.66, "MG": 1.41, "AL": 1.21, "K": 2.03, "CA": 1.76, "MN": 1.39,
    "FE": 1.32, "CO": 1.26, "NI": 1.24, "CU": 1.32, "ZN": 1.22, "MO": 1.54, "CD": 1.44,
    "HG": 1.32, "LI": 1.28, "PT": 1.36, "AU": 1.36, "AG": 1.45, "PD": 1.39,
}
_DEFAULT_RADIUS = 0.77


def covalent_radius(element: str) -> float:
    return COVALENT_RADII.get(element.upper(), _DEFAULT_RADIUS)


def build_graph(atoms, bond_tol: float = 0.45, skip_elements=("H",)) -> nx.Graph:
    """Heavy-atom graph of a list of :class:`posemap.io.Atom`.

    Nodes are indices into ``atoms`` (after skipping ``skip_elements``); each carries an
    ``el`` attribute (upper-case element). Edges connect atoms whose distance is below the
    summed covalent radii plus ``bond_tol`` (Å).
    """
    idx = [i for i, a in enumerate(atoms) if a.element.upper() not in {e.upper() for e in skip_elements}]
    G = nx.Graph()
    for i in idx:
        G.add_node(i, el=atoms[i].element.upper())
    xyz = {i: atoms[i].xyz for i in idx}
    rad = {i: covalent_radius(atoms[i].element) for i in idx}
    for a_pos, i in enumerate(idx):
        for j in idx[a_pos + 1:]:
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if d < rad[i] + rad[j] + bond_tol:
                G.add_edge(i, j)
    return G


def element_multiset(graph: nx.Graph) -> "dict[str, int]":
    counts: dict[str, int] = {}
    for _, el in graph.nodes(data="el"):
        counts[el] = counts.get(el, 0) + 1
    return counts
