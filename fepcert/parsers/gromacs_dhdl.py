"""
Parser for GROMACS dhdl.xvg files across multiple lambda windows.
"""

from typing import Dict, Any, List, Optional, Tuple
import os
import glob
import re
import numpy as np


def parse_gromacs_dhdl(filepath: str) -> Dict[str, Any]:
    """
    Parses a single GROMACS dhdl.xvg file.

    Parameters
    ----------
    filepath : str
        Path to dhdl.xvg file.

    Returns
    -------
    data : dict
        time, dH/dlambda gradients, and lambda state value.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    time_pts = []
    dhdl_pts = []
    lambda_val = None
    
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            if line_str.startswith("@"):
                # Regex match for lambda value: fep-lambda 0.2000 or vdw-lambda 0.2000
                m = re.search(r"(?:fep|vdw|coul|lambda)[-_]lambda\s+([0-9\.]+)", line_str, re.IGNORECASE)
                if m:
                    try:
                        lambda_val = float(m.group(1))
                    except ValueError:
                        pass
                continue
                
            parts = line_str.split()
            if len(parts) >= 2:
                try:
                    t = float(parts[0])
                    # column 1 is dH/dlambda in standard GROMACS dhdl.xvg
                    grad = float(parts[1])
                    time_pts.append(t)
                    dhdl_pts.append(grad)
                except ValueError:
                    pass

    return {
        "time": np.asarray(time_pts, dtype=float),
        "gradients": np.asarray(dhdl_pts, dtype=float),
        "lambda": lambda_val,
        "n_frames": len(dhdl_pts)
    }


def parse_gromacs_dhdl_directory(dirpath: str, file_pattern: str = "dhdl*.xvg") -> Dict[str, Any]:
    """
    Parses a collection of GROMACS dhdl.xvg files from a directory, sorted by lambda.

    Parameters
    ----------
    dirpath : str
        Directory path containing lambda window files.
    file_pattern : str, default 'dhdl*.xvg'

    Returns
    -------
    data : dict
        lambda_values (list of float), gradients_list (list of np.ndarray).
    """
    files = sorted(glob.glob(os.path.join(dirpath, file_pattern)))
    if not files:
        raise FileNotFoundError(f"No files matching pattern '{file_pattern}' found in {dirpath}")
        
    parsed_items = []
    for idx, f in enumerate(files):
        d = parse_gromacs_dhdl(f)
        lam = d["lambda"] if d["lambda"] is not None else float(idx / max(1, len(files) - 1))
        parsed_items.append((lam, d["gradients"]))

    # Sort by lambda
    parsed_items.sort(key=lambda x: x[0])
    
    lambdas = [item[0] for item in parsed_items]
    grads = [item[1] for item in parsed_items]

    return {
        "lambda_values": lambdas,
        "gradients_list": grads,
        "n_windows": len(lambdas)
    }
