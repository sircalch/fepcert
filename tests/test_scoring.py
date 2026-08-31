"""
Tests for scoring, reporting, and CLI execution in FEPCert.
"""

import os
import tempfile
import numpy as np
import pytest
from fepcert.core.scoring import assess_fep_quality
from fepcert.reporters.plot_generator import generate_fepcert_figures
from fepcert.reporters.manuscript_prep import generate_fepcert_manuscript_assets
from fepcert.reporters.html_report import generate_fepcert_html_report
from fepcert.cli import run_demo


def test_full_fep_validation_pipeline():
    meta = {
        "transformation": "LigA -> LigB",
        "engine": "GROMACS",
        "temperature_k": 298.15
    }
    
    lambdas = [0.0, 0.25, 0.5, 0.75, 1.0]
    rng = np.random.default_rng(42)
    # Good phase space overlap between windows
    grads = [rng.normal(-4.0 * l, 1.2, 300) for l in lambdas]
    
    report = assess_fep_quality(
        metadata=meta,
        lambda_values=lambdas,
        gradients_list=grads,
        unit="kcal/mol"
    )
    
    assert report.overall_status in ["PASS", "WARNING"]
    assert report.free_energy.method == "Thermodynamic Integration (TI)"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Vector figures
        plots = generate_fepcert_figures(report, tmpdir, formats=["png", "svg"])
        assert len(plots) > 0
        for p in plots:
            assert os.path.exists(p)
            
        # Manuscript assets
        assets = generate_fepcert_manuscript_assets(report, tmpdir)
        assert os.path.exists(assets["summary_csv"])
        assert os.path.exists(assets["summary_tex"])
        assert os.path.exists(assets["methods_text"])
        assert os.path.exists(assets["citation_bib"])
        
        # HTML report
        html_p = os.path.join(tmpdir, "report.html")
        generate_fepcert_html_report(report, html_p, methods_text="Sample methods", citation_bib="@software{}")
        assert os.path.exists(html_p)
        assert os.path.getsize(html_p) > 500


def test_cli_demo_execution():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_demo(output_dir=tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "report.html"))
        assert os.path.exists(os.path.join(tmpdir, "fepcert_summary_table.csv"))
        assert os.path.exists(os.path.join(tmpdir, "citation.bib"))
