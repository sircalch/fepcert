"""
Tests for time convergence and hysteresis evaluation.
"""

import numpy as np
import pytest
from fepcert.core.convergence import evaluate_time_convergence


def test_time_convergence_stationary():
    lambdas = [0.0, 0.5, 1.0]
    rng = np.random.default_rng(42)
    # Stationary series
    grads = [rng.normal(2.0 * l, 0.2, 500) for l in lambdas]
    
    res = evaluate_time_convergence(lambdas, grads, n_fractions=5)
    assert res.status == "PASS"
    assert res.final_hysteresis < 0.20
    assert res.dissipated_work < 0.10
