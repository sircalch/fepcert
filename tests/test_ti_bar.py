"""
Tests for Thermodynamic Integration (TI) and Bennett Acceptance Ratio (BAR).
"""

import numpy as np
import pytest
from fepcert.core.ti_bar import calculate_ti_free_energy, calculate_bar_free_energy


def test_ti_integration():
    lambdas = [0.0, 0.25, 0.5, 0.75, 1.0]
    # Linear function: y = 2 * lambda => integral_0^1 (2*x) dx = 1.0
    grads = [np.array([2.0 * l] * 100) for l in lambdas]
    
    res = calculate_ti_free_energy(lambdas, grads, unit="kcal/mol")
    assert np.isclose(res.delta_g, 1.0, atol=1e-3)
    assert res.status == "PASS"
    assert res.delta_g_error == 0.0


def test_bar_free_energy():
    lambdas = [0.0, 0.5, 1.0]
    # 2 windows
    # Window 0->1: Delta U around 2.0 kcal/mol
    w_fwd = [np.random.normal(2.0, 0.5, 500), np.random.normal(3.0, 0.5, 500)]
    w_rev = [np.random.normal(-2.0, 0.5, 500), np.random.normal(-3.0, 0.5, 500)]
    
    res = calculate_bar_free_energy(lambdas, w_fwd, w_rev, temperature_k=298.15)
    assert res.status == "PASS"
    assert np.isclose(res.delta_g, 5.0, atol=0.5)
