"""
Parser for chemical perturbation network CSV files (edges, DeltaDeltaG values, and errors).
"""

from typing import List, Dict, Any
import os
import pandas as pd


def parse_perturbation_network_csv(filepath: str) -> List[Dict[str, Any]]:
    """
    Parses a perturbation network CSV file.
    Expected columns: 'ligand_a', 'ligand_b', 'delta_g' (or 'ddG'), optional 'error'.

    Parameters
    ----------
    filepath : str
        Path to network CSV file.

    Returns
    -------
    edges : list of dict
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    df = pd.read_csv(filepath)
    # Standardize column names
    col_map = {}
    for c in df.columns:
        c_low = c.lower().strip()
        if c_low in ["ligand_a", "lig_a", "source", "from"]:
            col_map[c] = "ligand_a"
        elif c_low in ["ligand_b", "lig_b", "target", "to"]:
            col_map[c] = "ligand_b"
        elif c_low in ["delta_g", "ddg", "calc_ddg", "free_energy"]:
            col_map[c] = "delta_g"
        elif c_low in ["error", "std", "sem", "sigma"]:
            col_map[c] = "error"
            
    df = df.rename(columns=col_map)
    if "ligand_a" not in df.columns or "ligand_b" not in df.columns or "delta_g" not in df.columns:
        raise ValueError(f"Network CSV must contain columns 'ligand_a', 'ligand_b', and 'delta_g'. Found: {list(df.columns)}")

    edges = []
    for _, row in df.iterrows():
        edges.append({
            "ligand_a": str(row["ligand_a"]),
            "ligand_b": str(row["ligand_b"]),
            "delta_g": float(row["delta_g"]),
            "error": float(row.get("error", 0.15)) if "error" in row and pd.notnull(row["error"]) else 0.15
        })

    return edges
