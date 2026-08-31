"""
Publication-ready vector figures for alchemical free energy calculations.
"""

from typing import List, Optional
import os
import numpy as np
import matplotlib.pyplot as plt
from fepcert.core.scoring import FEPValidationReport

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'figure.dpi': 300,
    'lines.linewidth': 2.0,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})


def generate_fepcert_figures(
    report: FEPValidationReport,
    output_dir: str,
    formats: List[str] = ("png", "svg", "pdf")
) -> List[str]:
    """
    Generates publication figures: Overlap Matrix, TI Gradient Curve, Time Convergence.

    Parameters
    ----------
    report : FEPValidationReport
    output_dir : str
    formats : list of str

    Returns
    -------
    saved_paths : list of str
    """
    os.makedirs(output_dir, exist_ok=True)
    saved_files = []

    # 1. Phase Space Overlap Matrix Heatmap
    if report.overlap_result is not None:
        mat = report.overlap_result.overlap_matrix
        k = mat.shape[0]
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        
        im = ax.imshow(mat, cmap="Blues", vmin=0.0, vmax=1.0, origin="lower")
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r"Phase Space Overlap $\Pi_{ij}$", rotation=270, labelpad=15)
        
        # Annotate values if K <= 16
        if k <= 16:
            for i in range(k):
                for j in range(k):
                    val = mat[i, j]
                    color = "white" if val > 0.5 else "black"
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontsize=8)
                    
        ax.set_xticks(range(k))
        ax.set_yticks(range(k))
        ax.set_xlabel(r"$\lambda$-State $j$")
        ax.set_ylabel(r"$\lambda$-State $i$")
        ax.set_title(rf"Phase Space Overlap Matrix ({k} states) — Min Adjacent: {report.overlap_result.min_adjacent_overlap*100:.1f}%")
        
        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"fepcert_overlap_matrix.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    # 2. TI Gradient Curve <dH/dlambda>
    fe = report.free_energy
    if fe.mean_gradients is not None and fe.gradient_errors is not None:
        lams = np.asarray(fe.lambda_values)
        means = np.asarray(fe.mean_gradients)
        errs = np.asarray(fe.gradient_errors)
        
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(lams, means, "o-", color="#0284c7", label=r"$\langle \partial H / \partial \lambda \rangle$")
        ax.fill_between(lams, means - errs, means + errs, color="#0284c7", alpha=0.25, label=r"$\pm 1\sigma$ SEM")
        
        # Shaded integration area
        ax.fill_between(lams, 0, means, color="#38bdf8", alpha=0.10)
        
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel(r"Alchemical Coupling Parameter $\lambda$")
        ax.set_ylabel(r"$\langle \partial H / \partial \lambda \rangle$ (kcal/mol)")
        ax.set_title(rf"Thermodynamic Integration Curve — $\Delta G = {fe.delta_g:.2f} \pm {fe.delta_g_error:.2f}$ {fe.unit}")
        ax.grid(True)
        ax.legend(loc="best", frameon=True)
        
        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"fepcert_ti_gradient_curve.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    # 3. Time Convergence & Forward/Reverse Plot
    if report.convergence_result is not None:
        cr = report.convergence_result
        t_frac = np.asarray(cr.time_fractions) * 100.0
        fwd = np.asarray(cr.forward_delta_g)
        rev = np.asarray(cr.reverse_delta_g)
        
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(t_frac, fwd, "s-", color="#16a34a", label=r"Forward Cumulative $\Delta G(t)$")
        ax.plot(t_frac, rev, "^-", color="#dc2626", label=r"Reverse Cumulative $\Delta G(t)$")
        
        # Final value horizontal dashed line
        ax.axhline(fwd[-1], color="#16a34a", linestyle=":", alpha=0.7)
        
        ax.set_xlim(5, 105)
        ax.set_xlabel("Simulation Trajectory Evaluated (%)")
        ax.set_ylabel(r"Calculated $\Delta G$ (kcal/mol)")
        ax.set_title(rf"Time-Series Convergence & Hysteresis — $\Delta\Delta G_{{\mathrm{{hyst}}}} = {cr.final_hysteresis:.2f}$ kcal/mol")
        ax.grid(True)
        ax.legend(loc="best", frameon=True)
        
        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"fepcert_time_convergence.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    return saved_files
