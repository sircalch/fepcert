"""
Validation of FEPCert against alchemlyb on real GROMACS free-energy data (alchemtest).

Data: alchemtest.gmx.load_benzene(). Benzene decoupled from water in two legs (Coulomb, 5 windows;
VDW, 16 windows), 40 ns per window, 300 K, energies in kJ/mol.

For every leg:
  * TI: FEPCert parse_gromacs_dhdl_directory + calculate_ti_free_energy on the raw .xvg.bz2 files,
    compared with alchemlyb.estimators.TI on extract_dHdl (all samples, and after
    decorrelate_dhdl).
  * BAR: forward/reverse reduced work between neighbouring windows built from alchemlyb's
    extract_u_nk, passed to FEPCert calculate_bar_free_energy, and compared with
    alchemlyb.estimators.BAR on the same u_nk.

    pip install alchemlyb alchemtest
    python validation/validate_alchemtest.py [--out validation/results]
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import fepcert  # noqa: E402
from fepcert.core.ti_bar import calculate_bar_free_energy, calculate_ti_free_energy  # noqa: E402
from fepcert.parsers.gromacs_dhdl import parse_gromacs_dhdl_directory  # noqa: E402

import alchemlyb  # noqa: E402
from alchemlyb.estimators import BAR, TI  # noqa: E402
from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk  # noqa: E402
from alchemlyb.preprocessing import decorrelate_dhdl  # noqa: E402
from alchemtest.gmx import load_benzene  # noqa: E402

try:
    from loguru import logger
    logger.remove()
except Exception:  # pragma: no cover
    pass

R_KJ = 0.0083144626181532  # kJ/(mol K)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    T = 300.0
    kT = R_KJ * T
    root = os.path.dirname(os.path.dirname(os.path.dirname(load_benzene().data["VDW"][0])))
    rows = []
    for leg in ("Coulomb", "VDW"):
        d = parse_gromacs_dhdl_directory(os.path.join(root, leg), "dhdl.xvg.bz2", recursive=True)
        ti_fc = calculate_ti_free_energy(d["lambda_values"], d["gradients_list"], unit="kJ/mol")

        files = [d_file for d_file in d["files"]]
        dh = [extract_dHdl(f, T=T) for f in files]
        ti_al = TI().fit(pd.concat(dh).sort_index())
        ti_al_dec = TI().fit(pd.concat([decorrelate_dhdl(x) for x in dh]).sort_index())

        unk = [extract_u_nk(f, T=T) for f in files]
        bar_al = BAR().fit(pd.concat(unk).sort_index())
        # Neighbour work in kJ/mol from reduced potentials: state k -> k+1 sampled at k, and back.
        states = list(unk[0].columns)
        lam_of_file = [float(u.index.get_level_values(1)[0]) for u in unk]
        order = np.argsort(lam_of_file)
        unk_sorted = [unk[i] for i in order]
        col_sorted = sorted(states)
        wf, wr = [], []
        for k in range(len(unk_sorted) - 1):
            a, b = col_sorted[k], col_sorted[k + 1]
            wf.append(((unk_sorted[k][b] - unk_sorted[k][a]).to_numpy()) * kT)
            wr.append(((unk_sorted[k + 1][a] - unk_sorted[k + 1][b]).to_numpy()) * kT)
        bar_fc = calculate_bar_free_energy([float(c) for c in col_sorted], wf, wr, temperature_k=T,
                                           unit="kJ/mol")

        rows.append({
            "leg": leg, "n_windows": d["n_windows"], "component": d["component"],
            "unit_read": d["unit"], "temperature_read_K": d["temperature_k"],
            "TI_fepcert": ti_fc.delta_g, "TI_fepcert_se": ti_fc.delta_g_error,
            "TI_alchemlyb": float(ti_al.delta_f_.iloc[0, -1] * kT),
            "TI_alchemlyb_se": float(ti_al.d_delta_f_.iloc[0, -1] * kT),
            "TI_alchemlyb_decorrelated": float(ti_al_dec.delta_f_.iloc[0, -1] * kT),
            "TI_alchemlyb_decorrelated_se": float(ti_al_dec.d_delta_f_.iloc[0, -1] * kT),
            "BAR_fepcert": bar_fc.delta_g, "BAR_fepcert_se": bar_fc.delta_g_error,
            "BAR_alchemlyb": float(bar_al.delta_f_.iloc[0, -1] * kT),
            "BAR_alchemlyb_se": float(np.sqrt(np.sum(np.diag(bar_al.d_delta_f_.to_numpy(), 1) ** 2)) * kT),
        })

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out, "alchemtest_benzene.csv"), index=False)
    pd.set_option("display.width", 200)
    print(df.T.to_string())
    meta = {"fepcert": fepcert.__version__, "alchemlyb": alchemlyb.__version__,
            "dataset": "alchemtest.gmx.load_benzene", "temperature_K": T, "unit": "kJ/mol"}
    with open(os.path.join(args.out, "versions.json"), "w") as fh:
        json.dump(meta, fh, indent=2)


if __name__ == "__main__":
    main()
