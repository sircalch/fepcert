"""
Thermodynamic cycle detection and closure error analysis in chemical perturbation networks.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np


@dataclass
class CycleItem:
    cycle_nodes: List[str]  # e.g. ['Lig1', 'Lig2', 'Lig3']
    cycle_edges: List[Tuple[str, str]]
    closure_error: float  # kcal/mol
    status: str


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


def _find_simple_cycles_3_4(graph: Dict[str, Dict[str, float]]) -> List[List[str]]:
    """
    Finds simple 3-node and 4-node directed cycles in perturbation graph.
    """
    nodes = list(graph.keys())
    cycles = []
    seen = set()

    # 1. 3-node cycles (triangles)
    for u in nodes:
        for v in graph.get(u, {}):
            for w in graph.get(v, {}):
                if u in graph.get(w, {}):
                    c = [u, v, w]
                    c_canon = tuple(sorted(c))
                    if c_canon not in seen:
                        seen.add(c_canon)
                        cycles.append(c)

    # 2. 4-node cycles (quadrilaterals)
    for u in nodes:
        for v in graph.get(u, {}):
            for w in graph.get(v, {}):
                if w != u:
                    for x in graph.get(w, {}):
                        if x != v and u in graph.get(x, {}):
                            c = [u, v, w, x]
                            c_canon = tuple(sorted(c))
                            if c_canon not in seen:
                                seen.add(c_canon)
                                cycles.append(c)

    return cycles


def evaluate_cycle_closure(
    edges: List[Dict[str, Any]],
    max_pass_rmse: float = 0.50,
    max_warn_rmse: float = 0.80
) -> CycleClosureResult:
    """
    Evaluates thermodynamic cycle closure errors in a relative binding affinity network.

    Parameters
    ----------
    edges : list of dict
        List of transformation edges: {'ligand_a': 'L1', 'ligand_b': 'L2', 'delta_g': 1.25, 'error': 0.15}.
    max_pass_rmse : float, default 0.50 kcal/mol
    max_warn_rmse : float, default 0.80 kcal/mol

    Returns
    -------
    result : CycleClosureResult
    """
    # Build adjacency map with signed delta_g (L_a -> L_b: +dG, L_b -> L_a: -dG)
    graph: Dict[str, Dict[str, float]] = {}
    nodes_set = set()
    
    for e in edges:
        u = str(e["ligand_a"])
        v = str(e["ligand_b"])
        dg = float(e["delta_g"])
        nodes_set.add(u)
        nodes_set.add(v)
        
        if u not in graph:
            graph[u] = {}
        if v not in graph:
            graph[v] = {}
            
        graph[u][v] = dg
        graph[v][u] = -dg

    # Find closed cycles
    raw_cycles = _find_simple_cycles_3_4(graph)
    cycle_items: List[CycleItem] = []
    errors_list = []

    for c_nodes in raw_cycles:
        n_c = len(c_nodes)
        c_edges = []
        c_error = 0.0
        for i in range(n_c):
            u = c_nodes[i]
            v = c_nodes[(i + 1) % n_c]
            c_edges.append((u, v))
            c_error += graph[u][v]
            
        st = "PASS" if abs(c_error) <= max_pass_rmse else ("WARNING" if abs(c_error) <= max_warn_rmse else "FAIL")
        cycle_items.append(CycleItem(
            cycle_nodes=c_nodes,
            cycle_edges=c_edges,
            closure_error=float(c_error),
            status=st
        ))
        errors_list.append(c_error)

    n_cycles = len(cycle_items)
    if n_cycles > 0:
        rmse = float(np.sqrt(np.mean(np.array(errors_list)**2)))
        max_err = float(np.max(np.abs(errors_list)))
    else:
        rmse = 0.0
        max_err = 0.0

    if n_cycles == 0:
        status = "PASS"
        diag = "Tree-like perturbation network (no closed thermodynamic cycles detected)."
    elif rmse <= max_pass_rmse and max_err <= max_warn_rmse:
        status = "PASS"
        diag = f"Thermodynamically consistent perturbation network ({n_cycles} cycles evaluated, Cycle RMSE = {rmse:.2f} kcal/mol <= {max_pass_rmse:.2f})."
    elif rmse <= max_warn_rmse:
        status = "WARNING"
        diag = f"Moderate cycle closure strain ({n_cycles} cycles, Cycle RMSE = {rmse:.2f} kcal/mol, Max error = {max_err:.2f} kcal/mol)."
    else:
        status = "FAIL"
        diag = f"Severe thermodynamic cycle closure failure (Cycle RMSE = {rmse:.2f} kcal/mol > {max_warn_rmse:.2f}, Max error = {max_err:.2f} kcal/mol). Network contains inconsistent transformations."

    return CycleClosureResult(
        n_nodes=len(nodes_set),
        n_edges=len(edges),
        n_cycles=n_cycles,
        cycles=cycle_items,
        cycle_rmse=rmse,
        max_cycle_error=max_err,
        status=status,
        diagnostic_message=diag
    )
