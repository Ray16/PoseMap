"""Structure I/O: minimal, dependency-light parsers for PDB and mmCIF.

Only what an atom-mapping task needs is parsed: element, atom name, coordinates, and the
(model, chain, residue, hetero-flag) grouping used to isolate individual ligands. For heavy
production mmCIF work install ``gemmi`` and feed its atoms to :class:`Structure` directly;
the built-in mmCIF reader handles the common ``_atom_site`` loop but not multi-block edge cases.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np

# Common non-ligand residue names excluded from ligand auto-detection.
WATER = {"HOH", "DOD", "WAT", "H2O"}
COMMON_IONS = {
    "NA", "K", "CL", "MG", "CA", "ZN", "MN", "FE", "CU", "NI", "CO", "CD", "HG",
    "SO4", "PO4", "NO3", "ACT", "EDO", "GOL", "PEG", "DMS", "IOD", "BR", "FLC",
}


@dataclass
class Atom:
    """One parsed atom. ``element`` is upper-case; ``xyz`` is a length-3 float array."""

    element: str
    name: str
    xyz: np.ndarray
    chain: str = ""
    resname: str = ""
    resseq: str = ""
    hetatm: bool = False
    serial: int = 0
    model: int = 1
    altloc: str = ""

    @property
    def res_key(self) -> tuple:
        """Identity of the residue this atom belongs to (model, chain, resname, resseq)."""
        return (self.model, self.chain, self.resname, self.resseq)


@dataclass
class Structure:
    """A parsed structure: a flat list of :class:`Atom` plus convenience selectors."""

    atoms: list[Atom] = field(default_factory=list)
    source: str = ""

    def __len__(self) -> int:
        return len(self.atoms)

    def select(self, hetatm: Optional[bool] = None, chain: Optional[str] = None,
               resname: Optional[str] = None, resseq: Optional[str] = None,
               model: Optional[int] = None, heavy_only: bool = False) -> list[Atom]:
        out = []
        for a in self.atoms:
            if hetatm is not None and a.hetatm != hetatm:
                continue
            if chain is not None and a.chain != chain:
                continue
            if resname is not None and a.resname != resname:
                continue
            if resseq is not None and str(a.resseq) != str(resseq):
                continue
            if model is not None and a.model != model:
                continue
            if heavy_only and a.element == "H":
                continue
            out.append(a)
        return out

    def residues(self, atoms: Optional[Iterable[Atom]] = None) -> "dict[tuple, list[Atom]]":
        """Group atoms by residue key, preserving file order within each residue."""
        groups: dict[tuple, list[Atom]] = {}
        for a in (self.atoms if atoms is None else atoms):
            groups.setdefault(a.res_key, []).append(a)
        return groups

    def ligand_groups(self, exclude_water: bool = True, exclude_ions: bool = True,
                      min_heavy: int = 2, heavy_only: bool = True) -> "dict[tuple, list[Atom]]":
        """Candidate ligand residues: HETATM residues that are not water/ions and have at
        least ``min_heavy`` heavy atoms. Keyed by residue key; values are (heavy) atom lists."""
        out: dict[tuple, list[Atom]] = {}
        for key, ats in self.residues(self.select(hetatm=True)).items():
            resname = key[2]
            if exclude_water and resname in WATER:
                continue
            if exclude_ions and resname in COMMON_IONS:
                continue
            heavy = [a for a in ats if a.element != "H"]
            use = heavy if heavy_only else ats
            if len(heavy) < min_heavy:
                continue
            out[key] = use
        return out


# ---------------------------------------------------------------------------------------
# Element inference
# ---------------------------------------------------------------------------------------
_TWO_LETTER = {"CL", "BR", "NA", "MG", "AL", "SI", "CA", "MN", "FE", "CO", "NI",
               "CU", "ZN", "SE", "MO", "CD", "HG", "LI", "BE", "KR", "XE", "PT",
               "AU", "AG", "PD", "SN", "SB", "TE", "BA", "CS", "RB", "SR"}


def _element_from_name(name: str) -> str:
    s = "".join(c for c in name if c.isalpha()).upper()
    if len(s) >= 2 and s[:2] in _TWO_LETTER:
        return s[:2]
    return s[:1] if s else "X"


# ---------------------------------------------------------------------------------------
# PDB
# ---------------------------------------------------------------------------------------
def parse_pdb(text: str, source: str = "") -> Structure:
    atoms: list[Atom] = []
    model = 1
    for ln in text.splitlines():
        rec = ln[:6]
        if rec == "MODEL ":
            try:
                model = int(ln[10:14])
            except ValueError:
                model += 1
            continue
        if rec == "ENDMDL":
            model += 1
            continue
        if rec not in ("ATOM  ", "HETATM"):
            continue
        name = ln[12:16].strip()
        elem = ln[76:78].strip().upper() or _element_from_name(name)
        try:
            xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
        try:
            serial = int(ln[6:11])
        except ValueError:
            serial = len(atoms) + 1
        atoms.append(Atom(
            element=elem, name=name, xyz=xyz,
            chain=ln[21].strip(), resname=ln[17:20].strip(), resseq=ln[22:26].strip(),
            hetatm=(rec == "HETATM"), serial=serial, model=model, altloc=ln[16].strip(),
        ))
    return Structure(atoms=atoms, source=source)


# ---------------------------------------------------------------------------------------
# mmCIF (common _atom_site loop)
# ---------------------------------------------------------------------------------------
def parse_cif(text: str, source: str = "") -> Structure:
    lines = text.splitlines()
    i = 0
    atoms: list[Atom] = []
    n = len(lines)
    while i < n:
        ln = lines[i].strip()
        if ln == "loop_":
            # collect the tag block
            tags: list[str] = []
            j = i + 1
            while j < n and lines[j].strip().startswith("_"):
                tags.append(lines[j].strip())
                j += 1
            if not tags or not tags[0].startswith("_atom_site."):
                i = j
                continue
            cols = {t.split(".", 1)[1]: k for k, t in enumerate(tags)}
            # data rows until next loop_/#/tag/empty
            while j < n:
                row = lines[j]
                s = row.strip()
                if s == "" or s == "#" or s == "loop_" or s.startswith("_"):
                    break
                toks = _cif_split(row)
                if len(toks) < len(tags):
                    j += 1
                    continue

                def g(key, default=""):
                    k = cols.get(key)
                    return toks[k] if k is not None and k < len(toks) else default

                if g("group_PDB", "ATOM") not in ("ATOM", "HETATM"):
                    j += 1
                    continue
                try:
                    xyz = np.array([float(g("Cartn_x")), float(g("Cartn_y")), float(g("Cartn_z"))])
                except ValueError:
                    j += 1
                    continue
                name = g("label_atom_id") or g("auth_atom_id")
                name = name.strip('"')
                elem = (g("type_symbol") or _element_from_name(name)).upper()
                try:
                    model = int(g("pdbx_PDB_model_num", "1"))
                except ValueError:
                    model = 1
                try:
                    serial = int(g("id", str(len(atoms) + 1)))
                except ValueError:
                    serial = len(atoms) + 1
                atoms.append(Atom(
                    element=elem, name=name, xyz=xyz,
                    chain=(g("auth_asym_id") or g("label_asym_id")).strip(),
                    resname=(g("auth_comp_id") or g("label_comp_id")).strip(),
                    resseq=(g("auth_seq_id") or g("label_seq_id")).strip(),
                    hetatm=(g("group_PDB", "ATOM") == "HETATM"),
                    serial=serial, model=model, altloc=g("label_alt_id", "").strip(),
                ))
                j += 1
            i = j
            continue
        i += 1
    return Structure(atoms=atoms, source=source)


def _cif_split(row: str) -> list[str]:
    """Split a CIF data row honoring single/double-quoted tokens."""
    out: list[str] = []
    i, n = 0, len(row)
    while i < n:
        c = row[i]
        if c in " \t":
            i += 1
            continue
        if c in "'\"":
            q = c
            i += 1
            start = i
            while i < n and row[i] != q:
                i += 1
            out.append(row[start:i])
            i += 1
        else:
            start = i
            while i < n and row[i] not in " \t":
                i += 1
            out.append(row[start:i])
    return out


def load_structure(path: str) -> Structure:
    """Load a .pdb/.ent or .cif/.mmcif file (by extension; falls back to PDB)."""
    with open(path) as fh:
        text = fh.read()
    lower = path.lower()
    if lower.endswith((".cif", ".mmcif")):
        return parse_cif(text, source=path)
    return parse_pdb(text, source=path)
