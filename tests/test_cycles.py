"""
Tests for thermodynamic cycle closure analysis.
"""

import numpy as np
import pytest
from fepcert.core.cycles import evaluate_cycle_closure


def test_cycle_closure_consistent():
    edges = [
        {"ligand_a": "L1", "ligand_b": "L2", "delta_g": 2.00, "error": 0.1},
        {"ligand_a": "L2", "ligand_b": "L3", "delta_g": 3.00, "error": 0.1},
        {"ligand_a": "L3", "ligand_b": "L1", "delta_g": -5.05, "error": 0.1} # Error = 2 + 3 - 5.05 = -0.05
    ]
    res = evaluate_cycle_closure(edges)
    assert res.status == "PASS"
    assert res.n_cycles == 1
    assert np.isclose(res.cycle_rmse, 0.05, atol=1e-2)


def test_cycle_closure_inconsistent():
    edges = [
        {"ligand_a": "L1", "ligand_b": "L2", "delta_g": 2.00, "error": 0.1},
        {"ligand_a": "L2", "ligand_b": "L3", "delta_g": 3.00, "error": 0.1},
        {"ligand_a": "L3", "ligand_b": "L1", "delta_g": -2.00, "error": 0.1} # Error = 2 + 3 - 2 = 3.0 kcal/mol
    ]
    res = evaluate_cycle_closure(edges)
    assert res.status == "FAIL"
    assert res.cycle_rmse > 1.5
