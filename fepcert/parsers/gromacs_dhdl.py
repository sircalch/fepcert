"""
Parser for GROMACS dhdl.xvg files across multiple lambda windows.

GROMACS writes one dhdl.xvg per lambda state. Its header identifies the state, the
temperature and the meaning of every column, for example::

    @ subtitle "T = 300 (K) \\xl\\f{} state 3: fep-lambda = 0.2000"
    @ subtitle "T = 298 (K) \\xl\\f{} state 7: (coul-lambda, vdw-lambda) = (1.0000, 0.3000)"
    @ s0 legend "dH/d\\xl\\f{} fep-lambda = 0.2000"
    @ s1 legend "dH/d\\xl\\f{} coul-lambda = 1.0000"
    @ s2 legend "\\xD\\f{}H \\xl\\f{} to 0.0000"

dH/dlambda columns are therefore identified by their legends, not by position (a total-energy
column may precede them), and the lambda value is read from the header instead of being
inferred from file names.
"""

from typing import Dict, Any, List, Optional, Tuple
import bz2
import glob
import gzip
import lzma
import os
import re
import numpy as np

_NUM = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
_TEMP_RE = re.compile(r"T\s*=\s*(" + _NUM + r")\s*\(K\)")
_SINGLE_STATE_RE = re.compile(r"([A-Za-z]+)-lambda\s*=?\s*(" + _NUM + r")")
_MULTI_STATE_RE = re.compile(r"\(([^)]*-lambda[^)]*)\)\s*=\s*\(([^)]*)\)")
_LEGEND_RE = re.compile(r"^@\s*s(\d+)\s+legend\s+\"(.*)\"")
_DHDL_LEGEND_RE = re.compile(r"dH/d.*?([A-Za-z]+)-lambda(?:\s*=?\s*(" + _NUM + r"))?")
_UNIT_RE = re.compile(r"\((kJ/mol|kcal/mol)")


def _open_text(filepath: str):
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8", errors="ignore")
    if filepath.endswith(".bz2"):
        return bz2.open(filepath, "rt", encoding="utf-8", errors="ignore")
    if filepath.endswith(".xz"):
        return lzma.open(filepath, "rt", encoding="utf-8", errors="ignore")
    return open(filepath, "r", encoding="utf-8", errors="ignore")


def parse_gromacs_dhdl(filepath: str) -> Dict[str, Any]:
    """
    Parses a single GROMACS dhdl.xvg file (optionally .gz / .bz2 / .xz compressed).

    Parameters
    ----------
    filepath : str
        Path to dhdl.xvg file.

    Returns
    -------
    data : dict
        'time', 'components' ({name: {'lambda': float or None, 'gradients': ndarray}}),
        'temperature_k', 'unit' (as written by GROMACS, normally kJ/mol), 'n_frames', and for
        single-component files the convenience keys 'lambda' and 'gradients'.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    temperature = None
    unit = None
    state_lambdas: Dict[str, float] = {}
    dhdl_columns: List[Tuple[int, str, Optional[float]]] = []
    rows: List[List[float]] = []

    with _open_text(filepath) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("@"):
                if "subtitle" in s:
                    m = _TEMP_RE.search(s)
                    if m:
                        temperature = float(m.group(1))
                    mm = _MULTI_STATE_RE.search(s)
                    if mm:
                        names = [n.strip().replace("-lambda", "") for n in mm.group(1).split(",")]
                        vals = [float(v) for v in mm.group(2).split(",")]
                        state_lambdas.update(dict(zip(names, vals)))
                    else:
                        ms = _SINGLE_STATE_RE.search(s)
                        if ms:
                            state_lambdas[ms.group(1).lower()] = float(ms.group(2))
                elif "yaxis" in s and "label" in s:
                    mu = _UNIT_RE.search(s)
                    if mu:
                        unit = mu.group(1)
                else:
                    ml = _LEGEND_RE.match(s)
                    if ml:
                        legend = ml.group(2)
                        md = _DHDL_LEGEND_RE.search(legend)
                        if md:
                            comp = md.group(1).lower()
                            lam = float(md.group(2)) if md.group(2) is not None else None
                            dhdl_columns.append((int(ml.group(1)) + 1, comp, lam))
                continue
            parts = s.split()
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue

    if not rows:
        raise ValueError(f"No numeric data found in {filepath}")
    width = min(len(r) for r in rows)
    data = np.asarray([r[:width] for r in rows], dtype=float)

    if not dhdl_columns:
        # Minimal files without legends: time + one dH/dlambda column.
        dhdl_columns = [(1, "fep", None)]

    components: Dict[str, Dict[str, Any]] = {}
    for col, comp, lam in dhdl_columns:
        if col >= width:
            continue
        lam_value = state_lambdas.get(comp, lam)
        components[comp] = {"lambda": lam_value, "gradients": data[:, col].copy()}

    out: Dict[str, Any] = {
        "time": data[:, 0].copy(),
        "components": components,
        "temperature_k": temperature,
        "unit": unit or "kJ/mol",
        "n_frames": int(data.shape[0]),
    }
    if len(components) == 1:
        only = next(iter(components.values()))
        out["lambda"] = only["lambda"]
        out["gradients"] = only["gradients"]
    else:
        out["lambda"] = None
        out["gradients"] = None
    return out


def parse_gromacs_dhdl_directory(dirpath: str, file_pattern: str = "dhdl*.xvg*",
                                 recursive: bool = False) -> Dict[str, Any]:
    """
    Parses a set of GROMACS dhdl.xvg files (one per lambda window), sorted by the lambda
    value read from each file header.

    Parameters
    ----------
    dirpath : str
        Directory containing the lambda-window files.
    file_pattern : str, default 'dhdl*.xvg*'
        Glob pattern (compressed files are accepted).
    recursive : bool, default False
        Also search sub-directories (e.g. one directory per window: 0000/dhdl.xvg).

    Returns
    -------
    data : dict
        'lambda_values', 'gradients_list', 'n_windows', 'component', 'unit', 'temperature_k',
        'files'.

    Raises
    ------
    ValueError
        If a window's lambda value cannot be read from its header, if lambda values repeat,
        or if more than one lambda component changes across the windows (analyse each leg,
        e.g. Coulomb and van der Waals, separately).
    """
    pattern = os.path.join(dirpath, "**", file_pattern) if recursive else os.path.join(dirpath, file_pattern)
    files = sorted(glob.glob(pattern, recursive=recursive))
    if not files:
        raise FileNotFoundError(f"No files matching pattern '{file_pattern}' found in {dirpath}")

    parsed = [(f, parse_gromacs_dhdl(f)) for f in files]
    comp_names = set.intersection(*(set(p["components"]) for _, p in parsed))
    if not comp_names:
        raise ValueError("The dhdl files do not share a common dH/dlambda component.")

    missing = [f for f, p in parsed for c in comp_names if p["components"][c]["lambda"] is None]
    if missing:
        raise ValueError(
            "Could not read the lambda value from the header of: "
            + ", ".join(os.path.basename(m) for m in missing[:5])
            + ". Lambda is never inferred from file names; provide lambda values explicitly.")

    varying = sorted(c for c in comp_names
                     if len({round(p["components"][c]["lambda"], 10) for _, p in parsed}) > 1)
    if len(varying) > 1:
        raise ValueError(
            f"Lambda components {varying} change across these windows. Analyse each leg "
            f"(e.g. Coulomb and van der Waals) in a separate directory.")
    comp = varying[0] if varying else sorted(comp_names)[0]

    items = sorted(((p["components"][comp]["lambda"], p["components"][comp]["gradients"], f, p)
                    for f, p in parsed), key=lambda t: t[0])
    lambdas = [it[0] for it in items]
    if len(set(np.round(lambdas, 10))) != len(lambdas):
        raise ValueError(f"Repeated lambda values for component '{comp}': {lambdas}")

    units = {it[3]["unit"] for it in items}
    temps = {it[3]["temperature_k"] for it in items if it[3]["temperature_k"] is not None}
    return {
        "lambda_values": lambdas,
        "gradients_list": [it[1] for it in items],
        "n_windows": len(lambdas),
        "component": comp,
        "unit": units.pop() if len(units) == 1 else "kJ/mol",
        "temperature_k": temps.pop() if len(temps) == 1 else None,
        "files": [it[2] for it in items],
    }
