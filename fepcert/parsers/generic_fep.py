"""
Generic CSV/TSV table parsers for alchemical gradients and energy distributions.
"""

from typing import Dict, Any, List
import os
import pandas as pd
import numpy as np


def parse_generic_fep_csv(filepath: str) -> Dict[str, Any]:
    """
    Parses a CSV or TSV file containing lambda windows as columns or stacked table.
    Expects columns: 'lambda_0.0', 'lambda_0.1', ... or columns 'lambda', 'gradient'.

    Parameters
    ----------
    filepath : str
        Path to CSV/TSV file.

    Returns
    -------
    data : dict
        lambda_values (list of float), gradients_list (list of np.ndarray).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    sep = "\t" if filepath.endswith(".tsv") else ","
    df = pd.read_csv(filepath, sep=sep)
    
    # Check if columns are lambda values (e.g. '0.0', '0.1', '0.2' or 'lambda_0.0')
    lambda_cols = []
    for c in df.columns:
        c_clean = str(c).lower().replace("lambda_", "").replace("l_", "")
        try:
            val = float(c_clean)
            lambda_cols.append((val, c))
        except ValueError:
            pass

    if lambda_cols:
        lambda_cols.sort(key=lambda x: x[0])
        lambdas = [x[0] for x in lambda_cols]
        grads = [df[x[1]].dropna().values for x in lambda_cols]
        return {
            "lambda_values": lambdas,
            "gradients_list": grads,
            "n_windows": len(lambdas)
        }

    # Stacked format with 'lambda' and 'gradient'
    if "lambda" in df.columns and "gradient" in df.columns:
        grouped = df.groupby("lambda")
        lambdas = sorted(grouped.groups.keys())
        grads = [grouped.get_group(l)["gradient"].dropna().values for l in lambdas]
        return {
            "lambda_values": [float(l) for l in lambdas],
            "gradients_list": grads,
            "n_windows": len(lambdas)
        }

    raise ValueError(f"Unable to detect lambda columns in CSV: {filepath}. Expected column names with lambda numbers (e.g. '0.0', '0.1', '0.2') or 'lambda' and 'gradient' columns.")
