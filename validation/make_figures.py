"""
Figures and tables for the FEPCert v2 manuscript, built only from validation/results/.

    python validation/make_figures.py
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import halfnorm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
from fepcert.core.cycles import evaluate_cycle_closure  # noqa: E402
from benchmark_cycles import load_edges  # noqa: E402

RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
TAB = os.path.join(HERE, "tables")
MM = 1 / 25.4
DOUBLE = 174 * MM
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
NEW, OLD, REF, AUX = "#2a78d6", "#e87ba4", "#1baf7a", "#eda100"
VAR = {"full": ("o", NEW), "short": ("s", AUX), "sparse": ("^", REF), "very_sparse": ("D", OLD)}
VARNAME = {"full": "all windows", "short": "first 5% of samples", "sparse": "every 2nd window", "very_sparse": "every 4th window"}


def setup():
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 1.2, "lines.markersize": 4,
        "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42})


def panel(ax, letter, x=-0.17):
    ax.text(x, 1.03, f"({letter})", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", color=INK)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def fig_overlap(d):
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 68 * MM))
    lim = (1.5e-3, 1.0)
    for v, (mk, col) in VAR.items():
        x = d[d.variant == v]
        a.scatter(x.ov_mbar_min, x.ov_fepcert_min, s=14, marker=mk, color=col, lw=0, label=f"1.2, {VARNAME[v]}", zorder=3)
        y = x[x.ov_v110_min.notna()]
        a.scatter(y.ov_mbar_min, y.ov_v110_min, s=14, marker=mk, facecolor="none", edgecolor=INK2, lw=0.6, zorder=2)
    a.scatter([], [], s=14, marker="o", facecolor="none", edgecolor=INK2, lw=0.6, label="1.1.0 (dH/d$\\lambda$ histograms)")
    a.plot(lim, lim, color=INK2, lw=0.6, ls="--")
    a.axvline(0.03, color=OLD, lw=0.6, ls=":")
    a.axhline(0.03, color=OLD, lw=0.6, ls=":")
    a.set(xscale="log", yscale="log", xlim=lim, ylim=lim, xlabel="smallest neighbouring overlap, multi-state MBAR",
          ylabel="smallest neighbouring overlap, FEPCert")
    a.legend(loc="lower right", fontsize=6)
    panel(a, "a")

    ref = d[d.variant == "full"].set_index(["system", "leg"]).bar_alchemlyb
    for v, (mk, col) in VAR.items():
        x = d[d.variant == v]
        dev = np.abs(x.ti_fepcert.to_numpy() - ref.loc[list(zip(x.system, x.leg))].to_numpy())
        b.scatter(x.ov_mbar_min, dev, s=14, marker=mk, color=col, lw=0, label=VARNAME[v], zorder=3)
    b.axvline(0.03, color=OLD, lw=0.6, ls=":")
    b.set(xscale="log", xlim=lim, ylim=(-0.2, 7.0), xlabel="smallest neighbouring overlap, multi-state MBAR",
          ylabel="|TI $-$ BAR (all windows)| (kJ mol$^{-1}$)")
    b.legend(loc="upper right")
    panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig2_overlap")


def fig_cycles(c):
    zs = []
    for p in sorted(glob.glob(os.path.join(HERE, "data", "fepplus_edges", "*", "*.csv"))):
        zs += [x.z_score for x in evaluate_cycle_closure(load_edges(p)).cycles]
    z = np.abs(np.array(zs))
    z = z[np.isfinite(z)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 66 * MM), gridspec_kw={"width_ratios": [1, 1.15]})
    mx = max(c.ref_cycles.max(), 1) * 1.05
    a.scatter(c.ref_cycles, c.v110_cycles, s=12, marker="^", color=OLD, lw=0, label="1.1.0", zorder=2)
    a.scatter(c.ref_cycles, c.fc_cycles, s=12, marker="o", color=NEW, lw=0, label="1.2", zorder=3)
    a.plot([0, mx], [0, mx], color=INK2, lw=0.6, ls="--")
    a.set(xlim=(0, mx), ylim=(0, mx), xlabel="cycles of 3 and 4 edges (networkx)", ylabel="cycles found by FEPCert")
    a.legend(loc="upper left")
    panel(a, "a")

    xs = np.sort(z)
    b.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=NEW, lw=1.4, label=f"FEP+ cycles (n = {len(z)})")
    grid = np.linspace(0, 8, 300)
    b.plot(grid, halfnorm.cdf(grid), color=INK2, lw=1.0, ls="--", label="expected if edge errors were correct")
    s = np.median(z) / halfnorm.median()
    b.plot(grid, halfnorm.cdf(grid / s), color=AUX, lw=1.0, ls=":", label=f"edge errors $\\times$ {s:.1f}")
    b.axvline(1.96, color=OLD, lw=0.6, ls=":")
    b.set(xlim=(0, 8), ylim=(0, 1.0), xlabel="|closure error| / combined Bennett error", ylabel="fraction of cycles")
    b.legend(loc="lower right")
    panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig3_cycles")
    return z, s


def tables(d, c, z, s):
    kj = r" (kJ\,mol$^{-1}$)"
    rows = [(r"TI $\Delta G$" + kj, (d.ti_fepcert - d.ti_alchemlyb).abs().max(), "alchemlyb TI", len(d)),
            (r"BAR $\Delta G$" + kj, (d.bar_fepcert - d.bar_alchemlyb).abs().max(), "alchemlyb BAR", len(d)),
            ("BAR standard error" + kj, (d.bar_fepcert_se - d.bar_alchemlyb_se).abs().max(), "alchemlyb BAR", len(d)),
            ("neighbouring-state overlap", d.ov_dev_pymbar_pair.max(), "pymbar two-state MBAR", len(d)),
            (r"forward TI, first 10\%" + kj, (d.conv_fepcert_fwd10 - d.conv_alchemlyb_fwd10).abs().max(),
             r"alchemlyb forward\_backward\_convergence", int(d.conv_alchemlyb_fwd10.notna().sum()))]
    lines = [r"\begin{tabular}{lrrl}", r"\toprule", r"Quantity & analyses & max.\ abs.\ difference & reference \\", r"\midrule"]
    for name, v, ref, n in rows:
        lines.append(f"{name} & {n} & {v:.1e} & {ref} " + r"\\")
    lines.append(r"cycle closure RMSE (kcal\,mol$^{-1}$) & 92 maps & "
                 + f"{(c.fc_cycle_rmse - c.ref_cycle_rmse).abs().max():.1e} & networkx cycles " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_agreement.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    lines = [r"\begin{tabular}{lllrrrrll}", r"\toprule",
             r"System & leg & variant & windows & TI & BAR & min.\ overlap (MBAR) & overlap 1.1.0 / 1.2 & convergence 1.1.0 / 1.2 \\",
             r"\midrule"]
    for r in d.itertuples():
        o_old = r.ov_v110_status if isinstance(r.ov_v110_status, str) else "n/a"
        c_old = r.conv_v110_status if isinstance(r.conv_v110_status, str) else "n/a"
        c_new = r.conv_fepcert_status if isinstance(r.conv_fepcert_status, str) else "n/a"
        vname = VARNAME[r.variant].replace("%", r"\%")
        lines.append(f"{r.system.replace('_', ' ')} & {r.leg} & {vname} & {int(r.n_windows)} & {r.ti_fepcert:.2f} & "
                     f"{r.bar_fepcert:.2f} & {r.ov_mbar_min:.3f} & {o_old} / {r.ov_fepcert_status} & {c_old} / {c_new} & {r.tibar_status} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_s_analyses.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(RES, "cycle_z_summary.txt"), "w") as fh:
        fh.write(f"cycles with finite z: {len(z)}\nmedian |z|: {np.median(z):.3f}\nfraction |z|>1.96: {(z > 1.96).mean():.3f}\n"
                 f"fraction |z|>3: {(z > 3).mean():.3f}\nscale (median/half-normal median): {s:.3f}\n"
                 f"scale (95th percentile/1.96): {np.quantile(z, 0.95) / 1.96:.3f}\n")


def main():
    for p in (FIG, TAB):
        os.makedirs(p, exist_ok=True)
    setup()
    d = pd.read_csv(os.path.join(RES, "benchmark.csv"))
    c = pd.read_csv(os.path.join(RES, "cycles.csv"))
    fig_overlap(d)
    z, s = fig_cycles(c)
    tables(d, c, z, s)
    print("figures written")


if __name__ == "__main__":
    main()
