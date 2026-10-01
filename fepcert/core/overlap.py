"""
Phase-space overlap between neighbouring lambda states.

The overlap is the element O_ij of the MBAR overlap matrix O = W^T W N (Klimovich, Shirts & Mobley,
J. Comput.-Aided Mol. Des. 2015, 29, 397), the quantity for which the usual 0.03 threshold was
proposed. For a pair of neighbouring states it is computed exactly from the reduced energy
differences sampled in each state (forward and reverse work), using the two-state MBAR (= BAR)
weights. Histograms of dH/dlambda do not measure phase-space overlap: two windows can have similar
dH/dlambda distributions and no overlap in configuration space, or the reverse. When only dH/dlambda
is available the check is reported as NOT_APPLICABLE.
"""

from typing import List, Optional
from dataclasses import dataclass
import numpy as np
from scipy import optimize

KB = {"kcal/mol": 0.001987204259, "kj/mol": 0.008314462618}


@dataclass
class OverlapAnalysisResult:
    n_states: int
    overlap_matrix: np.ndarray  # tridiagonal: diagonal and neighbouring-pair elements only
    adjacent_overlaps: List[float]
    min_adjacent_overlap: float
    mean_adjacent_overlap: float
    status: str  # 'PASS', 'WARNING', 'FAIL', 'NOT_APPLICABLE'
    diagnostic_message: str


def _calculate_1d_distribution_overlap(samples_i: np.ndarray, samples_j: np.ndarray, n_bins: int = 50) -> float:
    r"""
    Histogram overlap integral \int min(p_i, p_j) of two one-dimensional samples. Kept for
    backward compatibility; it is not a measure of phase-space overlap between lambda states.
    """
    s_i = samples_i[np.isfinite(samples_i)]
    s_j = samples_j[np.isfinite(samples_j)]
    if len(s_i) < 5 or len(s_j) < 5:
        return 0.0
    u_min = min(np.min(s_i), np.min(s_j))
    u_max = max(np.max(s_i), np.max(s_j))
    if np.isclose(u_min, u_max):
        return 1.0
    bins = np.linspace(u_min, u_max, n_bins + 1)
    hist_i, _ = np.histogram(s_i, bins=bins, density=True)
    hist_j, _ = np.histogram(s_j, bins=bins, density=True)
    return float(np.clip(np.sum(np.minimum(hist_i, hist_j) * np.diff(bins)), 0.0, 1.0))


def two_state_overlap(w_forward_reduced: np.ndarray, w_reverse_reduced: np.ndarray) -> float:
    """
    MBAR overlap O_01 of two states from reduced work values (dimensionless, in units of kT):
    w_forward = u_1 - u_0 evaluated on samples of state 0, w_reverse = u_0 - u_1 on samples of state 1.
    """
    wf = np.asarray(w_forward_reduced, dtype=float)
    wr = np.asarray(w_reverse_reduced, dtype=float)
    wf, wr = wf[np.isfinite(wf)], wr[np.isfinite(wr)]
    n0, n1 = len(wf), len(wr)
    if n0 == 0 or n1 == 0:
        return 0.0
    a = np.concatenate([wf, -wr])  # u_1 - u_0 for every sample

    def log_denominator(df):  # log(N0 + N1 exp(df - a)), computed stably
        return np.logaddexp(np.log(n0), np.log(n1) + df - a)

    def g(df):  # sum_n W_n0 - 1 = 0 defines df (two-state MBAR / BAR)
        return np.sum(np.exp(-log_denominator(df))) - 1.0

    lo, hi = float(np.min(a)) - 50.0, float(np.max(a)) + 50.0
    df = optimize.brentq(g, lo, hi)
    log_w0 = -log_denominator(df)
    log_w1 = df - a + log_w0
    return float(n1 * np.sum(np.exp(log_w0 + log_w1)))


def calculate_overlap_matrix(
    energy_differences_list: Optional[List[np.ndarray]] = None,
    min_pass_overlap: float = 0.08,
    min_warn_overlap: float = 0.03,
    forward_work_list: Optional[List[np.ndarray]] = None,
    reverse_work_list: Optional[List[np.ndarray]] = None,
    temperature_k: float = 298.15,
    unit: str = "kcal/mol"
) -> OverlapAnalysisResult:
    """
    Overlap between neighbouring lambda states.

    Parameters
    ----------
    energy_differences_list : list of np.ndarray, optional
        dH/dlambda samples per window. Accepted for backward compatibility only: they do not
        determine phase-space overlap, and on their own the result is NOT_APPLICABLE.
    forward_work_list, reverse_work_list : list of np.ndarray
        For each neighbouring pair k, k+1: U_{k+1} - U_k on samples of state k, and U_k - U_{k+1}
        on samples of state k+1, in `unit`.
    temperature_k : float
    unit : 'kcal/mol' or 'kJ/mol'
    min_pass_overlap, min_warn_overlap : float
        Thresholds on the smallest neighbouring-pair overlap (0.03 after Klimovich et al. 2015).
    """
    if forward_work_list is None or reverse_work_list is None:
        k = len(energy_differences_list) if energy_differences_list is not None else 0
        return OverlapAnalysisResult(
            n_states=k, overlap_matrix=np.full((k, k), np.nan), adjacent_overlaps=[], min_adjacent_overlap=float("nan"),
            mean_adjacent_overlap=float("nan"), status="NOT_APPLICABLE",
            diagnostic_message=("Phase-space overlap needs the energy differences between neighbouring states "
                                "(forward and reverse work, or u_nk); dH/dlambda alone does not determine it."))
    beta = 1.0 / (KB[unit.lower()] * temperature_k)
    n_pairs = len(forward_work_list)
    k = n_pairs + 1
    adj = [two_state_overlap(np.asarray(f) * beta, np.asarray(r) * beta)
           for f, r in zip(forward_work_list, reverse_work_list)]
    o = np.eye(k)
    for i, v in enumerate(adj):
        o[i, i + 1] = o[i + 1, i] = v
    min_adj, mean_adj = float(np.min(adj)), float(np.mean(adj))
    bottlenecks = [f"{i}->{i + 1} ({v:.3f})" for i, v in enumerate(adj) if v < min_warn_overlap]
    if min_adj >= min_pass_overlap:
        status = "PASS"
        diag = f"Smallest neighbouring-state overlap {min_adj:.3f} (>= {min_pass_overlap:g})."
    elif min_adj >= min_warn_overlap:
        status = "WARNING"
        diag = f"Smallest neighbouring-state overlap {min_adj:.3f}; adding intermediate states would improve precision."
    else:
        status = "FAIL"
        diag = (f"Smallest neighbouring-state overlap {min_adj:.3f} < {min_warn_overlap:g}; free energies across "
                f"{', '.join(bottlenecks)} are unreliable. Add intermediate states.")
    return OverlapAnalysisResult(n_states=k, overlap_matrix=o, adjacent_overlaps=adj, min_adjacent_overlap=min_adj,
                                 mean_adjacent_overlap=mean_adj, status=status, diagnostic_message=diag)
