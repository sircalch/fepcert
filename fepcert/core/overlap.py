"""
Phase space overlap matrix calculation and Bennett/Klimovich criterion evaluation.
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np


@dataclass
class OverlapAnalysisResult:
    n_states: int
    overlap_matrix: np.ndarray
    adjacent_overlaps: List[float]
    min_adjacent_overlap: float
    mean_adjacent_overlap: float
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str


def _calculate_1d_distribution_overlap(samples_i: np.ndarray, samples_j: np.ndarray, n_bins: int = 50) -> float:
    r"""
    Computes histogram probability density overlap integral between two samples:
    \Pi_{ij} = \int \min(p_i(u), p_j(u)) du
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
    
    bin_widths = np.diff(bins)
    overlap = float(np.sum(np.minimum(hist_i, hist_j) * bin_widths))
    return float(np.clip(overlap, 0.0, 1.0))


def calculate_overlap_matrix(
    energy_differences_list: List[np.ndarray],
    min_pass_overlap: float = 0.08,
    min_warn_overlap: float = 0.03
) -> OverlapAnalysisResult:
    """
    Calculates phase space overlap matrix across all lambda intermediate states.

    Parameters
    ----------
    energy_differences_list : list of np.ndarray
        List of energy samples (e.g. forward Delta U or dH/dlambda) sampled at each lambda state k.
    min_pass_overlap : float, default 0.08
        Minimum required adjacent overlap Pi_{k, k+1} for robust convergence.
    min_warn_overlap : float, default 0.03
        Klimovich minimum overlap cutoff (3%).

    Returns
    -------
    result : OverlapAnalysisResult
        Calculated overlap matrix, adjacent overlaps, and certification status.
    """
    k = len(energy_differences_list)
    if k < 2:
        raise ValueError("At least 2 lambda states are required to calculate overlap matrix.")
        
    overlap_matrix = np.eye(k)
    adjacent_overlaps = []
    
    for i in range(k):
        for j in range(i + 1, k):
            ov = _calculate_1d_distribution_overlap(
                energy_differences_list[i],
                energy_differences_list[j]
            )
            overlap_matrix[i, j] = ov
            overlap_matrix[j, i] = ov
            if j == i + 1:
                adjacent_overlaps.append(ov)

    min_adj = float(np.min(adjacent_overlaps)) if adjacent_overlaps else 0.0
    mean_adj = float(np.mean(adjacent_overlaps)) if adjacent_overlaps else 0.0
    
    # Identify bottleneck pairs
    bottlenecks = []
    for idx, ov in enumerate(adjacent_overlaps):
        if ov < min_warn_overlap:
            bottlenecks.append(f"lambda {idx} -> {idx+1} (Pi = {ov*100:.1f}%)")

    if min_adj >= min_pass_overlap:
        status = "PASS"
        diag = f"Excellent phase space overlap across all lambda windows (Min adjacent Pi = {min_adj*100:.1f}% >= {min_pass_overlap*100:.1f}%)."
    elif min_adj >= min_warn_overlap:
        status = "WARNING"
        diag = f"Moderate phase space overlap (Min adjacent Pi = {min_adj*100:.1f}%). Insertion of intermediate lambda windows recommended."
    else:
        status = "FAIL"
        b_str = ", ".join(bottlenecks)
        diag = f"Insufficient phase space overlap (Min adjacent Pi = {min_adj*100:.1f}% < {min_warn_overlap*100:.1f}%). Bottlenecks: {b_str}. High risk of bias in BAR/MBAR."

    return OverlapAnalysisResult(
        n_states=k,
        overlap_matrix=overlap_matrix,
        adjacent_overlaps=adjacent_overlaps,
        min_adjacent_overlap=min_adj,
        mean_adjacent_overlap=mean_adj,
        status=status,
        diagnostic_message=diag
    )
