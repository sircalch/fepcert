"""
Thermodynamic cycle detection and closure error analysis in chemical perturbation networks.

Every simple cycle of the undirected perturbation graph up to `max_length` edges is enumerated,
each identified by its set of edges (not its set of nodes, which would merge the distinct cycles
that share four or more nodes). The closure error of a cycle is the sum of the signed edge free
energies around it. When edge uncertainties are given, the expected spread of the closure error is
sigma_c = sqrt(sum sigma_e^2) over the edges of the cycle, and z = closure / sigma_c tells whether
the closure is consistent with the reported statistical errors.
"""

from typing import List, Dict, Any, Tuple, Optional
from dataclasses import dataclass
import numpy as np


@dataclass
class CycleItem:
    cycle_nodes: List[str]  # e.g. ['Lig1', 'Lig2', 'Lig3']
    cycle_edges: List[Tuple[str, str]]
    closure_error: float  # kcal/mol
    status: str
    expected_sigma: float = float("nan")  # sqrt of the summed edge variances
    z_score: float = float("nan")  # closure_error / expected_sigma


@dataclass
class CycleClosureResult:
    n_nodes: int
    n_edges: int
    n_cycles: int
    cycles: List[CycleItem]
    cycle_rmse: float  # Root mean square cycle error
    max_cycle_error: float
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str
    n_duplicate_edges: int = 0
    frac_cycles_z_above_2: float = float("nan")  # fraction of cycles with |z| > 1.96


def _simple_cycles(adj: Dict[str, Dict[str, float]], max_length: int) -> List[List[str]]:
    """All simple cycles with 3..max_length nodes of an undirected graph, each listed once."""
    order = {n: i for i, n in enumerate(sorted(adj))}
    found, seen = [], set()

    def dfs(start, path, visited):
        u = path[-1]
        for v in adj[u]:
            if v == start and len(path) >= 3:
                key = frozenset(frozenset(e) for e in zip(path, path[1:] + [start]))
                if key not in seen:
                    seen.add(key)
                    found.append(list(path))
            elif v not in visited and order[v] > order[start] and len(path) < max_length:
                visited.add(v)
                path.append(v)
                dfs(start, path, visited)
                path.pop()
                visited.discard(v)

    for s in sorted(adj, key=order.get):  # each cycle is found from its lowest-ordered node
        dfs(s, [s], {s})
    return found


def evaluate_cycle_closure(
    edges: List[Dict[str, Any]],
    max_pass_rmse: float = 0.50,
    max_warn_rmse: float = 0.80,
    max_length: int = 4
) -> CycleClosureResult:
    """
    Evaluates thermodynamic cycle closure errors in a relative binding affinity network.

    Parameters
    ----------
    edges : list of dict
        {'ligand_a': 'L1', 'ligand_b': 'L2', 'delta_g': 1.25, 'error': 0.15}; 'error' is optional.
    max_pass_rmse, max_warn_rmse : float, kcal/mol
    max_length : int, default 4
        Longest cycle (number of edges) enumerated.

    Notes
    -----
    Repeated edges between the same pair of ligands are combined by inverse-variance weighting
    when errors are given (by a plain mean otherwise) and counted in `n_duplicate_edges`.
    """
    pairs: Dict[Tuple[str, str], List[Tuple[float, float]]] = {}
    for e in edges:
        u, v, dg = str(e["ligand_a"]), str(e["ligand_b"]), float(e["delta_g"])
        err = e.get("error")
        err = float(err) if err is not None and np.isfinite(float(err)) and float(err) > 0 else np.nan
        key = (u, v) if u < v else (v, u)
        pairs.setdefault(key, []).append((dg if key == (u, v) else -dg, err))

    adj: Dict[str, Dict[str, float]] = {}
    sig: Dict[Tuple[str, str], float] = {}
    n_dup = 0
    for (u, v), vals in pairs.items():
        n_dup += len(vals) - 1
        d = np.array([x for x, _ in vals])
        s = np.array([y for _, y in vals])
        if np.all(np.isfinite(s)):
            w = 1.0 / s ** 2
            dg, sg = float(np.sum(w * d) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w)))
        else:
            dg, sg = float(np.mean(d)), float("nan")
        adj.setdefault(u, {})[v] = dg
        adj.setdefault(v, {})[u] = -dg
        sig[(u, v)] = sig[(v, u)] = sg

    cycle_items: List[CycleItem] = []
    for nodes in _simple_cycles(adj, max_length):
        c_edges = [(nodes[i], nodes[(i + 1) % len(nodes)]) for i in range(len(nodes))]
        err = float(sum(adj[a][b] for a, b in c_edges))
        sc = float(np.sqrt(sum(sig[e] ** 2 for e in c_edges)))
        st = "PASS" if abs(err) <= max_pass_rmse else ("WARNING" if abs(err) <= max_warn_rmse else "FAIL")
        cycle_items.append(CycleItem(cycle_nodes=nodes, cycle_edges=c_edges, closure_error=err, status=st,
                                     expected_sigma=sc, z_score=err / sc if np.isfinite(sc) and sc > 0 else float("nan")))

    errors = np.array([c.closure_error for c in cycle_items])
    n_cycles = len(cycle_items)
    rmse = float(np.sqrt(np.mean(errors ** 2))) if n_cycles else 0.0
    max_err = float(np.max(np.abs(errors))) if n_cycles else 0.0
    z = np.array([c.z_score for c in cycle_items])
    frac_z = float(np.mean(np.abs(z[np.isfinite(z)]) > 1.96)) if np.isfinite(z).any() else float("nan")
    z_txt = f"; {frac_z:.0%} of cycles exceed 1.96 combined edge errors" if np.isfinite(frac_z) else ""

    if n_cycles == 0:
        status = "PASS"
        diag = "Tree-like perturbation network (no closed thermodynamic cycles detected)."
    elif rmse <= max_pass_rmse and max_err <= max_warn_rmse:
        status = "PASS"
        diag = f"{n_cycles} cycles, closure RMSE {rmse:.2f} kcal/mol (<= {max_pass_rmse:.2f}){z_txt}."
    elif rmse <= max_warn_rmse:
        status = "WARNING"
        diag = f"{n_cycles} cycles, closure RMSE {rmse:.2f} kcal/mol, largest {max_err:.2f} kcal/mol{z_txt}."
    else:
        status = "FAIL"
        diag = (f"{n_cycles} cycles, closure RMSE {rmse:.2f} kcal/mol > {max_warn_rmse:.2f}, largest {max_err:.2f} kcal/mol{z_txt}; "
                f"some transformations are inconsistent.")

    return CycleClosureResult(
        n_nodes=len(adj), n_edges=len(pairs), n_cycles=n_cycles, cycles=cycle_items, cycle_rmse=rmse,
        max_cycle_error=max_err, status=status, diagnostic_message=diag, n_duplicate_edges=n_dup,
        frac_cycles_z_above_2=frac_z)
