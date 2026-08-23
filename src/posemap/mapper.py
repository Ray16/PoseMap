"""Core: map a template molecule's atom identities onto a structure's (scrambled) atoms.

The mapping is a graph isomorphism between the template's heavy-atom graph (from a SMILES or
an RDKit mol, where every atom's chemical identity is known) and a ligand's heavy-atom graph
built from 3D coordinates. Because it matches *connectivity + element*, it is invariant to the
arbitrary atom names/order a docking or co-folding tool emits, and it fails cleanly (returns
``matched=False``) when a pose's connectivity does not match the template — which is exactly the
signal you want to gate out distorted poses.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher, categorical_node_match

from .io import Structure, Atom
from .graph import build_graph, element_multiset

try:
    from rdkit import Chem
except Exception:  # pragma: no cover - rdkit is a hard dep for template building
    Chem = None


@dataclass
class MappedAtom:
    """A pose atom identified with a template atom (by chemical identity)."""

    template_idx: int
    local_idx: int          # index into the ligand's heavy-atom list
    element: str
    name: str               # the (arbitrary) name it carries in the structure file
    xyz: np.ndarray
    serial: int

    def __repr__(self) -> str:
        return (f"MappedAtom(template_idx={self.template_idx}, name={self.name!r}, "
                f"element={self.element!r}, xyz=[{self.xyz[0]:.3f},{self.xyz[1]:.3f},{self.xyz[2]:.3f}])")


@dataclass
class MapResult:
    matched: bool
    reason: str = ""
    ligand_key: Optional[tuple] = None
    pose_atoms: list = field(default_factory=list)          # heavy Atom list of the ligand
    mappings: list = field(default_factory=list)            # list of {template_idx: local_idx}
    template: object = None                                  # RDKit mol
    candidates: list = field(default_factory=list)          # other ligand keys that also matched

    # -- queries ------------------------------------------------------------------------
    def atom(self, template_idx: int) -> list[MappedAtom]:
        """All pose atoms that any isomorphism maps template atom ``template_idx`` to.
        More than one only for symmetry-equivalent atoms."""
        locs = []
        seen = set()
        for mp in self.mappings:
            li = mp.get(template_idx)
            if li is not None and li not in seen:
                seen.add(li)
                locs.append(li)
        return [self._mapped(template_idx, li) for li in locs]

    def atoms_by_smarts(self, smarts: str, mapnum: int = 1) -> list[MappedAtom]:
        """Pose atoms corresponding to the template atom(s) bearing SMARTS map number
        ``mapnum``. Union over all SMARTS matches and all isomorphisms (symmetry)."""
        tidx = self.template_smarts_indices(smarts, mapnum)
        out: list[MappedAtom] = []
        seen = set()
        for ti in tidx:
            for ma in self.atom(ti):
                if ma.local_idx not in seen:
                    seen.add(ma.local_idx)
                    out.append(ma)
        return out

    def template_smarts_indices(self, smarts: str, mapnum: int = 1) -> set[int]:
        """Template atom indices carrying ``mapnum`` in any match of ``smarts``."""
        if self.template is None or Chem is None:
            return set()
        q = Chem.MolFromSmarts(smarts)
        if q is None:
            raise ValueError(f"bad SMARTS: {smarts}")
        qpos = [a.GetIdx() for a in q.GetAtoms() if a.GetAtomMapNum() == mapnum]
        if not qpos:
            raise ValueError(f"SMARTS has no atom with map number {mapnum}: {smarts}")
        qpos = qpos[0]
        return {m[qpos] for m in self.template.GetSubstructMatches(q)}

    def _mapped(self, template_idx: int, local_idx: int) -> MappedAtom:
        a: Atom = self.pose_atoms[local_idx]
        return MappedAtom(template_idx=template_idx, local_idx=local_idx,
                          element=a.element, name=a.name, xyz=a.xyz, serial=a.serial)

    def __bool__(self) -> bool:
        return self.matched


class PoseMapper:
    """Build once from a template; map onto many structures."""

    def __init__(self, template_graph: nx.Graph, template_mol=None, bond_tol: float = 0.45):
        self._SG = template_graph
        self.template = template_mol
        self.bond_tol = bond_tol
        self._multiset = element_multiset(template_graph)

    # -- constructors -------------------------------------------------------------------
    @classmethod
    def from_smiles(cls, smiles: str, bond_tol: float = 0.45) -> "PoseMapper":
        if Chem is None:
            raise RuntimeError("RDKit is required to build a template from SMILES")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"bad SMILES: {smiles}")
        return cls(_mol_graph(mol), template_mol=mol, bond_tol=bond_tol)

    @classmethod
    def from_mol(cls, mol, bond_tol: float = 0.45) -> "PoseMapper":
        """Template from an RDKit mol (e.g. from a CCD/mol block with correct bond orders)."""
        return cls(_mol_graph(mol), template_mol=mol, bond_tol=bond_tol)

    # -- mapping ------------------------------------------------------------------------
    def map(self, structure: Structure, chain: Optional[str] = None,
            resname: Optional[str] = None, resseq: Optional[str] = None,
            atoms: Optional[list] = None) -> MapResult:
        """Map the template onto one ligand of ``structure``.

        Specify the ligand with ``chain``/``resname``/``resseq`` or pass an explicit heavy-atom
        ``atoms`` list. With no selection, the ligand whose element composition matches the
        template is auto-detected (ambiguity is reported in ``result.candidates``).
        """
        if atoms is not None:
            groups = {("explicit",): [a for a in atoms if a.element != "H"]}
        elif chain is not None or resname is not None or resseq is not None:
            sel = structure.select(chain=chain, resname=resname, resseq=resseq, heavy_only=True)
            groups = structure.residues(sel) if (resseq is None and (chain or resname)) else {("sel",): sel}
            if not groups:
                return MapResult(False, reason="selection matched no atoms")
        else:
            groups = structure.ligand_groups()
            if not groups:
                return MapResult(False, reason="no ligand groups found for auto-detection")

        # rank candidate groups: element-composition match first
        matches_comp = [k for k, ats in groups.items()
                        if element_multiset(build_graph(ats, self.bond_tol)) == self._multiset]
        ordered = matches_comp + [k for k in groups if k not in matches_comp]

        best = None
        for key in ordered:
            ats = groups[key]
            res = self._map_one(ats, key)
            if res.matched:
                res.candidates = [k for k in matches_comp if k != key]
                return res
            if best is None:
                best = res
        # note: MapResult.__bool__ is `matched`, so use an explicit None check, not `or`
        return best if best is not None else MapResult(False, reason="no candidate ligand groups")

    def map_all(self, structure: Structure) -> list[MapResult]:
        """Map every ligand group whose composition matches the template (all copies)."""
        out = []
        for key, ats in structure.ligand_groups().items():
            if element_multiset(build_graph(ats, self.bond_tol)) != self._multiset:
                continue
            res = self._map_one(ats, key)
            if res.matched:
                out.append(res)
        return out

    def _map_one(self, atoms: list, key) -> MapResult:
        if len(atoms) != self._SG.number_of_nodes():
            return MapResult(False, reason=f"heavy-atom count {len(atoms)} != template "
                             f"{self._SG.number_of_nodes()}", ligand_key=key, pose_atoms=atoms,
                             template=self.template)
        PG = build_graph(atoms, self.bond_tol, skip_elements=())  # atoms already heavy
        # relabel pose graph nodes to 0..k local indices matching `atoms`
        gm = GraphMatcher(PG, self._SG, node_match=categorical_node_match("el", "C"))
        maps = []
        for mp in gm.isomorphisms_iter():                # mp: pose_local_idx -> template_idx
            maps.append({ti: li for li, ti in mp.items()})
        if not maps:
            return MapResult(False, reason="no graph isomorphism (connectivity mismatch)",
                             ligand_key=key, pose_atoms=atoms, template=self.template)
        return MapResult(True, ligand_key=key, pose_atoms=atoms, mappings=maps,
                         template=self.template)


def _mol_graph(mol) -> nx.Graph:
    G = nx.Graph()
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 1:
            continue
        G.add_node(a.GetIdx(), el=a.GetSymbol().upper())
    for b in mol.GetBonds():
        i, j = b.GetBeginAtom(), b.GetEndAtom()
        if i.GetAtomicNum() == 1 or j.GetAtomicNum() == 1:
            continue
        G.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    return G
