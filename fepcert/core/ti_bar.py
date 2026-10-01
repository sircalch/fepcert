"""
Thermodynamic Integration (TI) and Bennett Acceptance Ratio (BAR) free energy estimators.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np
from scipy import integrate, optimize


def statistical_inefficiency(series: np.ndarray, c_window: float = 6.0) -> float:
    """
    Statistical inefficiency g = 1 + 2 * sum_k C(k) of a correlated timeseries, using the
    Madras-Sokal self-consistent window (stop at the first negative C(k) or when
    k >= c_window * tau_int). Returns 1.0 for series shorter than 4 samples.
    """
    x = np.asarray(series, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 4:
        return 1.0
    xc = x - x.mean()
    var = float(np.dot(xc, xc) / n)
    if var == 0.0:
        return 1.0
    n_fft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
    f = np.fft.rfft(xc, n=n_fft)
    acov = np.fft.irfft(f * np.conjugate(f), n=n_fft)[:n] / (n - np.arange(n))
    acf = acov / var
    tau = 0.5
    for k in range(1, n // 2):
        if acf[k] < 0:
            break
        tau += acf[k]
        if k >= c_window * tau:
            break
    return max(1.0, 2.0 * tau)


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
            # Standard error of the mean corrected for autocorrelation: s * sqrt(g / n)
            sem = np.std(g_clean, ddof=1) * np.sqrt(statistical_inefficiency(g_clean) / n) if n > 1 else 0.0
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
    
    # Error propagation for the trapezoid rule. Delta G = sum_i w_i * <dH/dl>_i with
    # w_0 = dl_0/2, w_i = (dl_{i-1} + dl_i)/2 for interior windows, w_K = dl_{K-1}/2, so
    # var(Delta G) = sum_i w_i^2 * sigma_i^2 (interior windows enter both adjacent intervals).
    weights = np.zeros(k)
    weights[:-1] += d_lambda / 2.0
    weights[1:] += d_lambda / 2.0
    total_error = float(np.sqrt(np.sum((weights * errors) ** 2)))

    to_unit = 4.184 if unit.lower() == "kj/mol" else 1.0  # thresholds are defined in kcal/mol
    if total_error <= 0.30 * to_unit:
        status = "PASS"
        diag = f"TI integration converged with high precision (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit})."
    elif total_error <= 0.80 * to_unit:
        status = "WARNING"
        diag = f"TI integration has moderate statistical error (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit}). Extended sampling recommended."
    else:
        status = "FAIL"
        diag = f"Large TI integration uncertainty (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit} > {0.80 * to_unit:.2f}). Insufficient sampling."

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
            
        # Asymptotic BAR variance (Shirts, Bair, Hooker & Pande, PRL 2003, 91, 140601):
        # var(beta*dF) = (<f^2>_F/<f>_F^2 - 1)/n_F + (<f^2>_R/<f>_R^2 - 1)/n_R, f = Fermi function
        fermi_f = 1.0 / (1.0 + np.exp(w_fwd - (df_val * beta) + c_const))
        fermi_r = 1.0 / (1.0 + np.exp(w_rev + (df_val * beta) - c_const))
        mean_f, mean_r = np.mean(fermi_f), np.mean(fermi_r)

        var_f = (np.mean(fermi_f**2) / mean_f**2 - 1.0) / n_f if mean_f > 0 else np.inf
        var_r = (np.mean(fermi_r**2) / mean_r**2 - 1.0) / n_r if mean_r > 0 else np.inf

        var_total = (var_f + var_r) / (beta**2)
        err = np.sqrt(max(0.0, var_total))
        
        pair_errors.append(float(err))
        cumulative_dg.append(float(cumulative_dg[-1] + df_val))

    total_dg = float(cumulative_dg[-1])
    total_error = float(np.sqrt(np.sum(np.array(pair_errors)**2)))

    to_unit = 4.184 if unit.lower() == "kj/mol" else 1.0  # thresholds are defined in kcal/mol
    if total_error <= 0.25 * to_unit:
        status = "PASS"
        diag = f"BAR free energy calculation converged (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit})."
    elif total_error <= 0.60 * to_unit:
        status = "WARNING"
        diag = f"BAR free energy calculation has moderate uncertainty (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit})."
    else:
        status = "FAIL"
        diag = f"High BAR estimation error (Delta G = {total_dg:.2f} \u00b1 {total_error:.2f} {unit} > {0.60 * to_unit:.2f})."

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


def calculate_ti_free_energy_multi(
    lambda_matrix: List[List[float]],
    gradients_by_component: List[Dict[str, np.ndarray]],
    components: List[str],
    unit: str = "kcal/mol"
) -> FreeEnergyResult:
    """
    TI along a path in which several lambda components change (for example bonded restraints, then
    Coulomb, then van der Waals in one set of GROMACS windows).

    Delta G = sum over components c of the trapezoidal integral of <dH/dlambda_c> over lambda_c,
    with the windows in path order (increasing sum of the lambda components); this is what
    alchemlyb's TI estimator computes. Each component's standard error is corrected for
    autocorrelation, and the variances are added.

    Parameters
    ----------
    lambda_matrix : list of list of float
        One row per window, one column per component, in the order of `components`.
    gradients_by_component : list of dict
        One dict per window mapping each component to its dH/dlambda_c time series.
    components : list of str
    unit : str
    """
    lam = np.asarray(lambda_matrix, dtype=float)
    order = np.argsort(lam.sum(axis=1), kind="stable")
    lam = lam[order]
    grads = [gradients_by_component[i] for i in order]
    total, var = 0.0, 0.0
    cumulative = np.zeros(len(lam))
    for j, c in enumerate(components):
        means, errs = [], []
        for w in grads:
            g = np.asarray(w[c], dtype=float)
            g = g[np.isfinite(g)]
            means.append(float(np.mean(g)) if len(g) else 0.0)
            errs.append(float(np.std(g, ddof=1) * np.sqrt(statistical_inefficiency(g) / len(g))) if len(g) > 1 else 0.0)
        means, errs = np.asarray(means), np.asarray(errs)
        dl = np.diff(lam[:, j])
        steps = (means[:-1] + means[1:]) / 2.0 * dl
        cumulative[1:] += np.cumsum(steps)
        total += float(np.sum(steps))
        wts = np.zeros(len(lam))
        wts[:-1] += dl / 2.0
        wts[1:] += dl / 2.0
        var += float(np.sum((wts * errs) ** 2))
    err = float(np.sqrt(var))
    to_unit = 4.184 if unit.lower() == "kj/mol" else 1.0
    status = "PASS" if err <= 0.30 * to_unit else ("WARNING" if err <= 0.80 * to_unit else "FAIL")
    return FreeEnergyResult(
        method="Thermodynamic Integration (TI, multi-component)",
        delta_g=total, delta_g_error=err, unit=unit,
        lambda_values=[float(x) for x in lam.sum(axis=1)], cumulative_delta_g=cumulative.tolist(),
        mean_gradients=None, gradient_errors=None, status=status,
        diagnostic_message=f"TI over {len(components)} lambda components ({', '.join(components)}): "
                           f"Delta G = {total:.2f} \u00b1 {err:.2f} {unit}.")


def check_estimator_agreement(ti: FreeEnergyResult, bar: FreeEnergyResult, z_max: float = 3.0,
                              warn_kcal: float = 0.25, fail_kcal: float = 1.0) -> Dict[str, Any]:
    """
    Compares TI and BAR on the same windows. Trapezoidal TI carries a quadrature error when the
    windows are too sparse where <dH/dlambda> is curved, which neither its statistical error nor
    the overlap reveals; BAR does not have this error. A difference larger than z_max combined
    standard errors and than warn_kcal (fail_kcal) kcal/mol gives WARNING (FAIL).
    """
    to_unit = 4.184 if ti.unit.lower() == "kj/mol" else 1.0
    diff = ti.delta_g - bar.delta_g
    sigma = float(np.sqrt(ti.delta_g_error ** 2 + bar.delta_g_error ** 2))
    z = abs(diff) / sigma if sigma > 0 else float("inf")
    if z > z_max and abs(diff) > fail_kcal * to_unit:
        status = "FAIL"
    elif z > z_max and abs(diff) > warn_kcal * to_unit:
        status = "WARNING"
    else:
        status = "PASS"
    msg = (f"TI - BAR = {diff:.2f} {ti.unit} ({z:.1f} combined standard errors)"
           + ("; TI is affected by the lambda spacing (quadrature error): add windows where dH/dlambda is curved."
              if status != "PASS" else "."))
    return {"ti_minus_bar": float(diff), "z": float(z), "status": status, "diagnostic_message": msg}
