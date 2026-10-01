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


def test_thresholds_follow_the_unit():
    import numpy as np
    from fepcert.core.convergence import evaluate_time_convergence
    rng = np.random.default_rng(1)
    lam = [0.0, 0.5, 1.0]
    # first and second halves differ by about 0.3 kcal/mol in total
    g_kcal = [np.concatenate([rng.normal(0.0, 0.05, 200), rng.normal(0.3, 0.05, 200)]) for _ in lam]
    kcal = evaluate_time_convergence(lam, g_kcal, unit="kcal/mol")
    kj = evaluate_time_convergence(lam, [g * 4.184 for g in g_kcal], unit="kJ/mol")
    assert kcal.status == kj.status
