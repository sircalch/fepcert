"""
Time-series convergence, forward/reverse chunking analysis, and hysteresis audit.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np


@dataclass
class ConvergenceAnalysisResult:
    time_fractions: List[float]  # [0.1, 0.2, ..., 1.0]
    forward_delta_g: List[float]
    reverse_delta_g: List[float]
    final_hysteresis: float
    dissipated_work: float
    max_late_drift: float
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str


def evaluate_time_convergence(
    lambda_values: List[float],
    gradients_time_series: List[np.ndarray],
    n_fractions: int = 10,
    max_pass_hysteresis: float = 0.40,
    max_warn_hysteresis: float = 0.80
) -> ConvergenceAnalysisResult:
    """
    Evaluates cumulative forward and reverse TI convergence as a function of simulation time.

    Parameters
    ----------
    lambda_values : list of float
        Sorted lambda schedule.
    gradients_time_series : list of np.ndarray
        Time series arrays of dH/dlambda for each state.
    n_fractions : int, default 10
        Number of time evaluation checkpoints (10%, 20%, ..., 100%).
    max_pass_hysteresis : float, default 0.40 kcal/mol
    max_warn_hysteresis : float, default 0.80 kcal/mol

    Returns
    -------
    result : ConvergenceAnalysisResult
    """
    lambdas = np.asarray(lambda_values, dtype=float)
    d_lambda = np.diff(lambdas)
    
    fractions = np.linspace(0.1, 1.0, n_fractions)
    fwd_dgs = []
    rev_dgs = []
    
    # Minimum trajectory length across windows
    min_len = min(len(g) for g in gradients_time_series)
    
    for f in fractions:
        n_pts = max(5, int(f * min_len))
        
        # Forward chunk (0 to n_pts)
        fwd_means = [np.mean(g[:n_pts]) for g in gradients_time_series]
        fwd_mid = (np.array(fwd_means[:-1]) + np.array(fwd_means[1:])) / 2.0
        fwd_dg = float(np.sum(fwd_mid * d_lambda))
        fwd_dgs.append(fwd_dg)
        
        # Reverse chunk (min_len - n_pts to min_len)
        rev_means = [np.mean(g[-n_pts:]) for g in gradients_time_series]
        rev_mid = (np.array(rev_means[:-1]) + np.array(rev_means[1:])) / 2.0
        rev_dg = float(np.sum(rev_mid * d_lambda))
        rev_dgs.append(rev_dg)

    # Hysteresis between two disjoint halves of the data (first half vs second half).
    # At 100% time the forward and reverse chunks contain identical samples, so their
    # difference is zero by construction and carries no information.
    half = max(1, min_len // 2)
    first_means = [np.mean(g[:half]) for g in gradients_time_series]
    second_means = [np.mean(g[min_len - half:min_len]) for g in gradients_time_series]
    dg_first = float(np.sum((np.array(first_means[:-1]) + np.array(first_means[1:])) / 2.0 * d_lambda))
    dg_second = float(np.sum((np.array(second_means[:-1]) + np.array(second_means[1:])) / 2.0 * d_lambda))
    final_hysteresis = abs(dg_first - dg_second)
    dissipated_work = 0.5 * final_hysteresis
    
    # Max drift across final 40% of time
    late_idx = int(n_fractions * 0.6)
    late_vals = fwd_dgs[late_idx:]
    max_drift = float(np.max(late_vals) - np.min(late_vals)) if len(late_vals) > 1 else 0.0

    if final_hysteresis <= max_pass_hysteresis and max_drift <= max_pass_hysteresis:
        status = "PASS"
        diag = f"High time-convergence stability (Hysteresis = {final_hysteresis:.2f} kcal/mol, Late drift = {max_drift:.2f} kcal/mol)."
    elif final_hysteresis <= max_warn_hysteresis:
        status = "WARNING"
        diag = f"Moderate convergence drift (Hysteresis = {final_hysteresis:.2f} kcal/mol). Transformation is approaching plateau but has mild fluctuations."
    else:
        status = "FAIL"
        diag = f"Severe convergence hysteresis ({final_hysteresis:.2f} kcal/mol > {max_warn_hysteresis:.2f} kcal/mol). Trajectory is trapped in non-equilibrated conformational states."

    return ConvergenceAnalysisResult(
        time_fractions=fractions.tolist(),
        forward_delta_g=fwd_dgs,
        reverse_delta_g=rev_dgs,
        final_hysteresis=final_hysteresis,
        dissipated_work=dissipated_work,
        max_late_drift=max_drift,
        status=status,
        diagnostic_message=diag
    )
