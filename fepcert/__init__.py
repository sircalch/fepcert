"""
FEPCert: Automated Quality-Control, Phase Space Overlap, and Thermodynamic
Cycle Closure Certification for Alchemical Free Energy Simulations (FEP, TI, BAR, MBAR).
"""

__version__ = "1.0.0"
__author__ = "Andre Monreal-Hernández"
__license__ = "MIT"

from fepcert.core.overlap import (
    calculate_overlap_matrix,
    OverlapAnalysisResult
)
from fepcert.core.ti_bar import (
    calculate_ti_free_energy,
    calculate_bar_free_energy,
    FreeEnergyResult
)
from fepcert.core.convergence import (
    evaluate_time_convergence,
    ConvergenceAnalysisResult
)
from fepcert.core.cycles import (
    evaluate_cycle_closure,
    CycleClosureResult
)
from fepcert.core.scoring import assess_fep_quality, FEPValidationReport

__all__ = [
    "__version__",
    "calculate_overlap_matrix",
    "OverlapAnalysisResult",
    "calculate_ti_free_energy",
    "calculate_bar_free_energy",
    "FreeEnergyResult",
    "evaluate_time_convergence",
    "ConvergenceAnalysisResult",
    "evaluate_cycle_closure",
    "CycleClosureResult",
    "assess_fep_quality",
    "FEPValidationReport"
]
