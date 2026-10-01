"""
Cycle-closure check of FEPCert on the 92 perturbation maps of the public FEP+ benchmark
(Ross et al. 2023; github.com/schrodinger/public_binding_free_energy_benchmark, commit a65d86b,
21_4_results/edge_predictions, MIT licence), using the Bennett (BAR) edge free energies.

Reference: every simple cycle of length 3 and 4 of the undirected perturbation graph enumerated by
networkx (simple_cycles with length_bound), with the closure error summed along the cycle.

    python validation/benchmark_cycles.py [--legacy DIR]

--legacy points to a checkout of FEPCert 1.1.0 to evaluate it on the same maps.
Output: validation/results/cycles.csv
"""
import argparse
import glob
import importlib
import os
import sys

import networkx as nx
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "fepplus_edges")


def load_edges(path):
    d = pd.read_csv(path).rename(columns={"Ligand1": "Lig 1", "Ligand2": "Lig 2"})
    return [{"ligand_a": str(r["Lig 1"]), "ligand_b": str(r["Lig 2"]), "delta_g": float(r["Bennett ddG (kcal/mol)"]),
             "error": float(r["Bennett std. error (kcal/mol)"])} for _, r in d.iterrows()]


def reference(edges, max_len=4):
    g = nx.Graph()
    dup = 0
    for e in edges:
        if g.has_edge(e["ligand_a"], e["ligand_b"]):
            dup += 1
        g.add_edge(e["ligand_a"], e["ligand_b"], dg={(e["ligand_a"], e["ligand_b"]): e["delta_g"],
                                                     (e["ligand_b"], e["ligand_a"]): -e["delta_g"]})
    errs = []
    for cyc in nx.simple_cycles(g, length_bound=max_len):
        tot = 0.0
        for i in range(len(cyc)):
            u, v = cyc[i], cyc[(i + 1) % len(cyc)]
            tot += g[u][v]["dg"][(u, v)]
        errs.append(tot)
    errs = np.abs(errs)
    return {"ref_nodes": g.number_of_nodes(), "ref_edges": g.number_of_edges(), "ref_duplicate_edges": dup,
            "ref_cycles": len(errs), "ref_cycle_rmse": float(np.sqrt(np.mean(errs ** 2))) if len(errs) else 0.0,
            "ref_cycle_max": float(errs.max()) if len(errs) else 0.0,
            "ref_cycle_basis": len(nx.cycle_basis(g))}


def fepcert_cycles(module_root, edges, prefix):
    sys.path.insert(0, module_root)
    for m in [m for m in sys.modules if m.startswith("fepcert")]:
        del sys.modules[m]
    cyc = importlib.import_module("fepcert.core.cycles")
    sys.path.pop(0)
    r = cyc.evaluate_cycle_closure(edges)
    out = {f"{prefix}_cycles": r.n_cycles, f"{prefix}_cycle_rmse": r.cycle_rmse, f"{prefix}_cycle_max": r.max_cycle_error,
           f"{prefix}_status": r.status}
    if hasattr(r, "frac_cycles_z_above_2"):
        z = np.array([c.z_score for c in r.cycles])
        out.update({f"{prefix}_frac_z_above_2": r.frac_cycles_z_above_2,
                    f"{prefix}_median_abs_z": float(np.median(np.abs(z))) if len(z) else np.nan,
                    f"{prefix}_n_cycles_z": int(np.isfinite(z).sum())})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legacy", default=None)
    a = ap.parse_args()
    rows = []
    for path in sorted(glob.glob(os.path.join(DATA, "*", "*.csv"))):
        edges = load_edges(path)
        row = {"group": os.path.basename(os.path.dirname(path)), "map": os.path.basename(path)[:-4], "n_edges_file": len(edges)}
        row.update(reference(edges))
        row.update(fepcert_cycles(os.path.dirname(HERE), edges, "fc"))
        if a.legacy:
            row.update(fepcert_cycles(a.legacy, edges, "v110"))
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(HERE, "results", "cycles.csv"), index=False)
    print(len(d), "maps;", int(d.ref_cycles.sum()), "reference cycles")
    for p in ("fc", "v110"):
        if f"{p}_cycles" in d:
            print(p, "cycles found", int(d[f"{p}_cycles"].sum()), "; maps with the reference cycle count",
                  int((d[f"{p}_cycles"] == d.ref_cycles).sum()), "; max |rmse diff|",
                  float((d[f"{p}_cycle_rmse"] - d.ref_cycle_rmse).abs().max()))


if __name__ == "__main__":
    main()
