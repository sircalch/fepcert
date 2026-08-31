"""
Tests for GROMACS and CSV parsers in FEPCert.
"""

import os
import tempfile
import numpy as np
import pytest
from fepcert.parsers.gromacs_dhdl import parse_gromacs_dhdl
from fepcert.parsers.generic_fep import parse_generic_fep_csv
from fepcert.parsers.network_csv import parse_perturbation_network_csv


def test_gromacs_dhdl_parser():
    content = """# GROMACS dhdl output
@ title "dH/dlambda"
@ xaxis label "Time (ps)"
@ yaxis label "dH/dlambda (kJ/mol)"
@ s0 legend "dH/dlambda fep-lambda 0.2000"
   0.0000   -12.4500
   2.0000   -11.8000
   4.0000   -13.1000
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".xvg", delete=False) as f:
        f.write(content)
        f_path = f.name

    try:
        data = parse_gromacs_dhdl(f_path)
        assert data["n_frames"] == 3
        assert np.isclose(data["lambda"], 0.2000)
        assert len(data["gradients"]) == 3
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)


def test_generic_csv_parser():
    content = """lambda_0.0,lambda_0.5,lambda_1.0
-5.2,-3.1,1.4
-4.9,-2.8,1.8
-5.1,-3.0,1.5
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        f_path = f.name

    try:
        data = parse_generic_fep_csv(f_path)
        assert data["n_windows"] == 3
        assert data["lambda_values"] == [0.0, 0.5, 1.0]
        assert len(data["gradients_list"]) == 3
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)


def test_network_csv_parser():
    content = """ligand_a,ligand_b,delta_g,error
Lig1,Lig2,-2.5,0.15
Lig2,Lig3,1.2,0.10
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        f_path = f.name

    try:
        edges = parse_perturbation_network_csv(f_path)
        assert len(edges) == 2
        assert edges[0]["ligand_a"] == "Lig1"
        assert edges[0]["delta_g"] == -2.5
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)
