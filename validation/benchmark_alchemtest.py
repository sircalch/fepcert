"""
Benchmark of FEPCert against alchemlyb/pymbar on real GROMACS free-energy data (alchemtest).

    pip install alchemlyb alchemtest
    python validation/benchmark_alchemtest.py --legacy <checkout of FEPCert 1.1.0>

Legs (alchemtest.gmx): benzene and ethanol hydration (Coulomb, VDW; one lambda component each),
the ligand and complex legs of an absolute binding free energy of n-phenylglycinonitrile to T4
lysozyme, and a water particle between surfaces (two or three lambda components changing
along one set of windows). For each leg, four variants:
    full         all windows, all samples
    short        the first 5% of the samples of every window
    sparse       every second window removed (end states kept)
    very_sparse  every fourth window kept (end states kept)
References (alchemlyb 2.5, pymbar 4):
    TI          alchemlyb.estimators.TI on extract_dHdl (multi-component lambda handled by alchemlyb)
    BAR         alchemlyb.estimators.BAR on extract_u_nk
    overlap     MBAR overlap matrix of all retained states (alchemlyb MBAR.overlap_matrix), and the
                two-state MBAR overlap of each neighbouring pair (pymbar), which is the quantity
                FEPCert 1.2 computes
    convergence alchemlyb.convergence.forward_backward_convergence (TI), first 10% forward/backward
Output: validation/results/benchmark.csv
"""
import argparse
import importlib
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from fepcert.core.ti_bar import (calculate_ti_free_energy, calculate_bar_free_energy, calculate_ti_free_energy_multi,  # noqa: E402
                                check_estimator_agreement)
from fepcert.core.overlap import calculate_overlap_matrix  # noqa: E402
from fepcert.core.convergence import evaluate_time_convergence  # noqa: E402
from fepcert.parsers.gromacs_dhdl import parse_gromacs_dhdl  # noqa: E402

import alchemtest.gmx as agmx  # noqa: E402
import pymbar  # noqa: E402
from alchemlyb.estimators import BAR, MBAR, TI  # noqa: E402
from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk  # noqa: E402
from alchemlyb.preprocessing import decorrelate_dhdl  # noqa: E402
from alchemlyb.convergence import forward_backward_convergence  # noqa: E402

try:
    from loguru import logger
    logger.remove()
except Exception:  # pragma: no cover
    pass

R_KJ = 0.0083144626181532
LEGS = [("benzene", "load_benzene", "Coulomb"), ("benzene", "load_benzene", "VDW"),
        ("ethanol", "load_ethanol", "Coulomb"), ("ethanol", "load_ethanol", "VDW"),
        ("ABFE", "load_ABFE", "ligand"), ("ABFE", "load_ABFE", "complex"),
        ("water_particle", "load_water_particle_with_total_energy", "AllStates")]
VARIANTS = ("full", "short", "sparse", "very_sparse")
LEGACY = {}


def windows(files):
    """FEPCert parser on each file; windows in path order, varying components, lambda matrix."""
    parsed = [parse_gromacs_dhdl(f) for f in files]
    comps = sorted(set.intersection(*(set(p["components"]) for p in parsed)))
    varying = [c for c in comps if len({round(p["components"][c]["lambda"], 10) for p in parsed}) > 1]
    order = sorted(range(len(parsed)), key=lambda i: sum(parsed[i]["components"][c]["lambda"] for c in varying))
    lam = [[parsed[i]["components"][c]["lambda"] for c in varying] for i in order]
    grads = [{c: parsed[i]["components"][c]["gradients"] for c in varying} for i in order]
    return order, varying, lam, grads


def state_of(u):
    key = u.index[0][1:]
    return key[0] if len(key) == 1 else tuple(key)


def pymbar_pair_overlap(u0, u1, a, b):
    """Two-state MBAR overlap of states a, b from the reduced potentials of their samples."""
    u_kn = np.vstack([np.concatenate([u0[a].to_numpy(), u1[a].to_numpy()]),
                      np.concatenate([u0[b].to_numpy(), u1[b].to_numpy()])])
    m = pymbar.MBAR(u_kn, np.array([len(u0), len(u1)]))
    return float(m.compute_overlap()["matrix"][0, 1])


def one_variant(name, leg, variant, files, T):
    kT = R_KJ * T
    row = {"system": name, "leg": leg, "variant": variant}
    order, comps, lam, grads = windows(files)
    files = [files[i] for i in order]
    n = len(files)
    idx = list(range(n))
    if variant == "sparse":
        idx = sorted(set(range(0, n, 2)) | {n - 1})
    elif variant == "very_sparse":
        idx = sorted(set(range(0, n, 4)) | {n - 1})
    frac = 0.05 if variant == "short" else 1.0

    def cut(x):
        return x[: max(20, int(len(x) * frac))]

    lam = [lam[i] for i in idx]
    grads = [{c: cut(grads[i][c]) for c in comps} for i in idx]
    dh = [cut(extract_dHdl(files[i], T=T)) for i in idx]
    unk = [cut(extract_u_nk(files[i], T=T)) for i in idx]
    cols = [state_of(u) for u in unk]
    unk = [u[cols] for u in unk]
    row.update(components="+".join(comps), n_components=len(comps), n_windows=len(idx),
               n_samples=int(np.median([len(g[comps[0]]) for g in grads])))

    # --- TI
    if len(comps) == 1:
        c = comps[0]
        ti = calculate_ti_free_energy([x[0] for x in lam], [g[c] for g in grads], unit="kJ/mol")
    else:
        ti = calculate_ti_free_energy_multi(lam, grads, comps, unit="kJ/mol")
    ti_al = TI().fit(pd.concat(dh).sort_index())
    ti_dec = TI().fit(pd.concat([decorrelate_dhdl(x) for x in dh]).sort_index())
    row.update(ti_fepcert=ti.delta_g, ti_fepcert_se=ti.delta_g_error, ti_fepcert_status=ti.status,
               ti_alchemlyb=float(ti_al.delta_f_.iloc[0, -1] * kT),
               ti_alchemlyb_dec=float(ti_dec.delta_f_.iloc[0, -1] * kT),
               ti_alchemlyb_dec_se=float(ti_dec.d_delta_f_.iloc[0, -1] * kT))

    # --- BAR (FEPCert) from neighbour work, and alchemlyb BAR / MBAR on the same u_nk
    wf = [(unk[k][cols[k + 1]] - unk[k][cols[k]]).to_numpy() * kT for k in range(len(unk) - 1)]
    wr = [(unk[k + 1][cols[k]] - unk[k + 1][cols[k + 1]]).to_numpy() * kT for k in range(len(unk) - 1)]
    bar = calculate_bar_free_energy(list(range(len(unk))), wf, wr, temperature_k=T, unit="kJ/mol")
    u_all = pd.concat(unk).sort_index()
    bar_al = BAR().fit(u_all)
    row.update(bar_fepcert=bar.delta_g, bar_fepcert_se=bar.delta_g_error,
               bar_alchemlyb=float(bar_al.delta_f_.iloc[0, -1] * kT),
               bar_alchemlyb_se=float(np.sqrt(np.sum(np.diag(bar_al.d_delta_f_.to_numpy(), 1) ** 2)) * kT))
    agree = check_estimator_agreement(ti, bar)
    row.update(tibar_status=agree["status"], tibar_z=agree["z"])

    # --- overlap
    ov = calculate_overlap_matrix(forward_work_list=wf, reverse_work_list=wr, temperature_k=T, unit="kJ/mol")
    ref_pair = [pymbar_pair_overlap(unk[k], unk[k + 1], cols[k], cols[k + 1]) for k in range(len(unk) - 1)]
    row.update(ov_fepcert_min=ov.min_adjacent_overlap, ov_fepcert_status=ov.status,
               ov_dev_pymbar_pair=float(np.max(np.abs(np.array(ov.adjacent_overlaps) - np.array(ref_pair)))),
               ov_pair_min=float(np.min(ref_pair)))
    try:
        O = MBAR().fit(u_all).overlap_matrix
        adj = [float(O[k, k + 1]) for k in range(len(unk) - 1)]
        row.update(ov_mbar_min=float(np.min(adj)), ov_mbar_n_below_003=int(np.sum(np.array(adj) < 0.03)))
    except Exception as e:  # noqa: BLE001
        row["mbar_error"] = str(e)[:120]
    if "overlap" in LEGACY and len(comps) == 1:
        old = LEGACY["overlap"].calculate_overlap_matrix([g[comps[0]] for g in grads])
        row.update(ov_v110_min=float(old.min_adjacent_overlap), ov_v110_status=old.status,
                   ov_v110_rank_corr=float(pd.Series(old.adjacent_overlaps).corr(pd.Series(ref_pair), method="spearman"))
                   if len(ref_pair) > 2 else np.nan)

    # --- time convergence (single-component legs)
    if len(comps) == 1:
        g1 = [g[comps[0]] for g in grads]
        lam1 = [x[0] for x in lam]
        conv = evaluate_time_convergence(lam1, g1, unit="kJ/mol")
        row.update(conv_fepcert_hysteresis=float(conv.final_hysteresis), conv_fepcert_status=conv.status,
                   conv_fepcert_fwd10=float(conv.forward_delta_g[0]), conv_fepcert_bwd10=float(conv.reverse_delta_g[0]))
        if "convergence" in LEGACY:
            old = LEGACY["convergence"].evaluate_time_convergence(lam1, g1)
            row.update(conv_v110_status=old.status)
        try:
            fb = forward_backward_convergence(dh, "TI", num=10)
            row.update(conv_alchemlyb_fwd10=float(fb.Forward.iloc[0] * kT), conv_alchemlyb_bwd10=float(fb.Backward.iloc[0] * kT))
        except Exception as e:  # noqa: BLE001
            row["conv_alchemlyb_error"] = str(e)[:120]
    return row


def load_legacy(path):
    saved = {m: sys.modules.pop(m) for m in list(sys.modules) if m.startswith("fepcert")}
    sys.path.insert(0, path)
    LEGACY["overlap"] = importlib.import_module("fepcert.core.overlap")
    LEGACY["convergence"] = importlib.import_module("fepcert.core.convergence")
    for m in [m for m in sys.modules if m.startswith("fepcert")]:
        del sys.modules[m]
    sys.modules.update(saved)
    sys.path.pop(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legacy", default=None)
    a = ap.parse_args()
    if a.legacy:
        load_legacy(a.legacy)
    rows = []
    for name, loader, leg in LEGS:
        ds = getattr(agmx, loader)()
        files = list(ds.data[leg])
        T = float(parse_gromacs_dhdl(files[0])["temperature_k"] or 300.0)
        for variant in VARIANTS:
            try:
                r = one_variant(name, leg, variant, files, T)
            except Exception as e:  # noqa: BLE001
                r = {"system": name, "leg": leg, "variant": variant, "error": f"{type(e).__name__}: {e}"[:200]}
            r["temperature_K"] = T
            rows.append(r)
            print(name, leg, variant, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                                       if k.startswith(("ti_f", "ti_alchemlyb", "bar_f", "bar_al", "ov_", "conv_", "error", "n_w"))},
                  flush=True)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "results", "benchmark.csv"), index=False)


if __name__ == "__main__":
    main()
