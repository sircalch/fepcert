"""
Thermodynamic Integration (TI) and Bennett Acceptance Ratio (BAR) free energy estimators.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np
from scipy import integrate, optimize


@dataclass
class FreeEnergyResult:
    method: str  # 'TI', 'BAR', 'MBAR'
    delta_g: float  # kcal/mol or kJ/mol
    delta_g_error: float
    unit: str
    lambda_values: List[float]
    cumulative_delta_g: List[float]
    mean_gradients: Optional[List[float]]
    gradient_errors: Optional[List[float]]
    status: str
    diagnostic_message: str


def calculate_ti_free_energy(
    lambda_values: List[float],
    gradients_list: List[np.ndarray],
    unit: str = "kcal/mol"
) -> FreeEnergyResult:
    """
    Computes free energy difference using Thermodynamic Integration (TI) with trapezoidal quadrature.

    Parameters
    ----------
    lambda_values : list of float
        Sorted lambda schedule from 0.0 to 1.0.
    gradients_list : list of np.ndarray
        List of dH/dlambda time-series samples for each lambda window.
    unit : str, default 'kcal/mol'

    Returns
    -------
    result : FreeEnergyResult
        Total Delta G, error, and cumulative integration curve.
    """
    lambdas = np.asarray(lambda_values, dtype=float)
    k = len(lambdas)
    
    means = []
    errors = []
    
    for grads in gradients_list:
        g = np.asarray(grads, dtype=float)
        g_clean = g[np.isfinite(g)]
        n = len(g_clean)
        if n == 0:
            means.append(0.0)
            errors.append(0.0)
        else:
            m = np.mean(g_clean)
            # Standard error of the mean
            sem = np.std(g_clean, ddof=1) / np.sqrt(n) if n > 1 else 0.0
            means.append(float(m))
            errors.append(float(sem))

    means = np.asarray(means)
    errors = np.asarray(errors)
    
    # Trapezoidal cumulative integration
    d_lambda = np.diff(lambdas)
    mid_means = (means[:-1] + means[1:]) / 2.0
    d_g_steps = mid_means * d_lambda
    
    cumulative_dg = [0.0]
    for step in d_g_steps:
        cumulative_dg.append(float(cumulative_dg[-1] + step))

    total_dg = float(cumulative_dg[-1])
    
    # Error propagation for trapezoid rule: sigma^2 = sum ( Delta\lambda_i / 2 )^2 * (sigma_i^2 + sigma_{i+1}^2)
    var_terms = (d_lambda / 2.0)**2 * (errors[:-1]**2 + errors[1:]**2)
    total_error = float(np.sqrt(np.sum(var_terms)))

    if total_error <= 0.30:
        status = "PASS"
        diag = f"TI integration converged with high precision (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit})."
    elif total_error <= 0.80:
        status = "WARNING"
        diag = f"TI integration has moderate statistical error (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit}). Extended sampling recommended."
    else:
        status = "FAIL"
        diag = f"Large TI integration uncertainty (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit} > 0.80). Insufficient sampling."

    return FreeEnergyResult(
        method="Thermodynamic Integration (TI)",
        delta_g=total_dg,
        delta_g_error=total_error,
        unit=unit,
        lambda_values=lambdas.tolist(),
        cumulative_delta_g=cumulative_dg,
        mean_gradients=means.tolist(),
        gradient_errors=errors.tolist(),
        status=status,
        diagnostic_message=diag
    )


def calculate_bar_free_energy(
    lambda_values: List[float],
    forward_work_list: List[np.ndarray],
    reverse_work_list: List[np.ndarray],
    temperature_k: float = 298.15,
    unit: str = "kcal/mol"
) -> FreeEnergyResult:
    """
    Computes free energy differences between adjacent lambda states using Bennett Acceptance Ratio (BAR).

    Parameters
    ----------
    lambda_values : list of float
        Sorted lambda schedule.
    forward_work_list : list of np.ndarray
        Forward Delta U_{k -> k+1} sampled at state k.
    reverse_work_list : list of np.ndarray
        Reverse Delta U_{k+1 -> k} sampled at state k+1.
    temperature_k : float, default 298.15
    unit : str, default 'kcal/mol'

    Returns
    -------
    result : FreeEnergyResult
    """
    # kB in kcal/(mol*K) or kJ/(mol*K)
    kb = 0.001987204259 if unit.lower() == "kcal/mol" else 0.008314462618
    beta = 1.0 / (kb * temperature_k)
    
    n_pairs = len(forward_work_list)
    cumulative_dg = [0.0]
    pair_errors = []
    
    for idx in range(n_pairs):
        w_fwd = np.asarray(forward_work_list[idx], dtype=float) * beta
        w_rev = np.asarray(reverse_work_list[idx], dtype=float) * beta
        
        n_f = len(w_fwd)
        n_r = len(w_rev)
        
        if n_f == 0 or n_r == 0:
            cumulative_dg.append(cumulative_dg[-1])
            pair_errors.append(0.0)
            continue
            
        # Self-consistent BAR equation
        c_const = np.log(n_f / n_r)
        
        def bar_root(df):
            lhs = np.sum(1.0 / (1.0 + np.exp(w_fwd - df + c_const)))
            rhs = np.sum(1.0 / (1.0 + np.exp(w_rev + df - c_const)))
            return lhs - rhs
            
        try:
            sol = optimize.root_scalar(bar_root, bracket=[-100.0, 100.0], method="brentq")
            df_val = sol.root / beta
        except Exception:
            # Fallback to mean difference
            df_val = float(np.mean(w_fwd)) / beta
            
        # Asymptotic BAR variance
        denom_f = np.mean(1.0 / (1.0 + np.exp(w_fwd - (df_val * beta) + c_const)))
        denom_r = np.mean(1.0 / (1.0 + np.exp(w_rev + (df_val * beta) - c_const)))
        
        var_f = (1.0 / max(1e-6, denom_f) - 1.0) / n_f if denom_f > 0 else 0.0
        var_r = (1.0 / max(1e-6, denom_r) - 1.0) / n_r if denom_r > 0 else 0.0
        
        var_total = (var_f + var_r) / (beta**2)
        err = np.sqrt(max(0.0, var_total))
        
        pair_errors.append(float(err))
        cumulative_dg.append(float(cumulative_dg[-1] + df_val))

    total_dg = float(cumulative_dg[-1])
    total_error = float(np.sqrt(np.sum(np.array(pair_errors)**2)))

    if total_error <= 0.25:
        status = "PASS"
        diag = f"BAR free energy calculation converged (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit})."
    elif total_error <= 0.60:
        status = "WARNING"
        diag = f"BAR free energy calculation has moderate uncertainty (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit})."
    else:
        status = "FAIL"
        diag = f"High BAR estimation error (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit} > 0.60)."

    return FreeEnergyResult(
        method="Bennett Acceptance Ratio (BAR)",
        delta_g=total_dg,
        delta_g_error=total_error,
        unit=unit,
        lambda_values=lambda_values,
        cumulative_delta_g=cumulative_dg,
        mean_gradients=None,
        gradient_errors=None,
        status=status,
        diagnostic_message=diag
    )
