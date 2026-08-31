"""
Alchemical Free Energy Validation Matrix, Decision Engine, and Certification Report.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass, asdict
import numpy as np

from fepcert.core.overlap import OverlapAnalysisResult, calculate_overlap_matrix
from fepcert.core.ti_bar import FreeEnergyResult, calculate_ti_free_energy, calculate_bar_free_energy
from fepcert.core.convergence import ConvergenceAnalysisResult, evaluate_time_convergence
from fepcert.core.cycles import CycleClosureResult, evaluate_cycle_closure


@dataclass
class FEPValidationReport:
    overall_status: str  # 'PASS', 'WARNING', 'FAIL'
    validation_score: str
    metadata: Dict[str, Any]
    free_energy: FreeEnergyResult
    overlap_result: Optional[OverlapAnalysisResult]
    convergence_result: Optional[ConvergenceAnalysisResult]
    cycle_result: Optional[CycleClosureResult]
    recommendations: List[str]
    provenance: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def assess_fep_quality(
    metadata: Dict[str, Any],
    lambda_values: List[float],
    gradients_list: Optional[List[np.ndarray]] = None,
    forward_work_list: Optional[List[np.ndarray]] = None,
    reverse_work_list: Optional[List[np.ndarray]] = None,
    network_edges: Optional[List[Dict[str, Any]]] = None,
    temperature_k: float = 298.15,
    unit: str = "kcal/mol"
) -> FEPValidationReport:
    """
    Evaluates alchemical free energy simulation convergence, lambda overlap, and cycle closure.

    Parameters
    ----------
    metadata : dict
        System description, perturbation (e.g. L1 -> L2), simulation engine (GROMACS, AMBER, OpenMM).
    lambda_values : list of float
        Sorted lambda schedule.
    gradients_list : list of np.ndarray, optional
        dH/dlambda samples per window for TI.
    forward_work_list : list of np.ndarray, optional
        Forward work samples for BAR/MBAR.
    reverse_work_list : list of np.ndarray, optional
        Reverse work samples for BAR/MBAR.
    network_edges : list of dict, optional
        Edges in chemical perturbation network for cycle analysis.
    temperature_k : float, default 298.15
    unit : str, default 'kcal/mol'

    Returns
    -------
    report : FEPValidationReport
        Comprehensive certification report.
    """
    statuses = []
    recommendations = []

    # 1. Free Energy Calculation (TI or BAR)
    if gradients_list is not None and len(gradients_list) > 0:
        fe_res = calculate_ti_free_energy(lambda_values, gradients_list, unit=unit)
    elif forward_work_list is not None and reverse_work_list is not None:
        fe_res = calculate_bar_free_energy(lambda_values, forward_work_list, reverse_work_list, temperature_k=temperature_k, unit=unit)
    else:
        raise ValueError("Must provide either gradients_list for TI or forward/reverse work lists for BAR.")

    statuses.append(fe_res.status)
    if fe_res.status != "PASS":
        recommendations.append(fe_res.diagnostic_message)

    # 2. Phase Space Overlap Analysis
    overlap_res = None
    samples_for_overlap = gradients_list if gradients_list is not None else forward_work_list
    if samples_for_overlap is not None and len(samples_for_overlap) >= 2:
        overlap_res = calculate_overlap_matrix(samples_for_overlap)
        statuses.append(overlap_res.status)
        if overlap_res.status != "PASS":
            recommendations.append(overlap_res.diagnostic_message)

    # 3. Time-Series Convergence Analysis
    conv_res = None
    if gradients_list is not None and len(gradients_list) > 0 and len(gradients_list[0]) >= 20:
        conv_res = evaluate_time_convergence(lambda_values, gradients_list)
        statuses.append(conv_res.status)
        if conv_res.status != "PASS":
            recommendations.append(conv_res.diagnostic_message)

    # 4. Thermodynamic Cycle Closure
    cycle_res = None
    if network_edges is not None and len(network_edges) >= 3:
        cycle_res = evaluate_cycle_closure(network_edges)
        statuses.append(cycle_res.status)
        if cycle_res.status != "PASS":
            recommendations.append(cycle_res.diagnostic_message)

    # Overall Decision
    if "FAIL" in statuses:
        overall_status = "FAIL"
        validation_score = "ALCHEMICAL FREE ENERGY = UNCONVERGED / HIGH ERROR (REJECTED)"
    elif "WARNING" in statuses:
        overall_status = "WARNING"
        validation_score = "ALCHEMICAL FREE ENERGY = ACCEPTABLE WITH WARNINGS"
    else:
        overall_status = "PASS"
        validation_score = "ALCHEMICAL FREE ENERGY = FULLY CERTIFIED (HIGH CONFIDENCE)"

    return FEPValidationReport(
        overall_status=overall_status,
        validation_score=validation_score,
        metadata=metadata,
        free_energy=fe_res,
        overlap_result=overlap_res,
        convergence_result=conv_res,
        cycle_result=cycle_res,
        recommendations=recommendations,
        provenance={
            "tool": "FEPCert",
            "version": "1.0.0",
            "citation": "Monreal-Hernández, A. (2026). FEPCert: An Open-Source Toolkit for Quality-Control, Phase Space Overlap, and Thermodynamic Cycle Closure Certification of Alchemical Free Energy Simulations."
        }
    )
