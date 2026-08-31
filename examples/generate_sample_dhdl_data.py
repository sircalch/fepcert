"""
Generates synthetic GROMACS dhdl.xvg files and perturbation network CSV for testing and tutorials.
"""

import os
import numpy as np


def generate_sample_fep_data(output_dir: str = "sample_gromacs_dhdl"):
    os.makedirs(output_dir, exist_ok=True)
    lambdas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    rng = np.random.default_rng(42)
    n_frames = 500
    
    # 1. Write dhdl files
    for idx, lam in enumerate(lambdas):
        fname = os.path.join(output_dir, f"dhdl_{idx:02d}_lambda_{lam:.2f}.xvg")
        base_grad = -10.0 * (lam ** 1.5) + 5.0 * lam - 3.0
        
        with open(fname, "w", encoding="utf-8") as f:
            f.write("# GROMACS dhdl output\n")
            f.write(f"@ s0 legend \"dH/dlambda fep-lambda {lam:.4f}\"\n")
            for t in range(n_frames):
                time_ps = t * 2.0
                grad = base_grad + rng.normal(0, 1.0)
                f.write(f"  {time_ps:10.4f}  {grad:12.6f}\n")
                
    print(f"Generated 11 sample dhdl.xvg files in: {os.path.abspath(output_dir)}/")

    # 2. Write network CSV
    net_file = os.path.join(output_dir, "sample_perturbation_network.csv")
    with open(net_file, "w", encoding="utf-8") as f:
        f.write("ligand_a,ligand_b,delta_g,error\n")
        f.write("Lig1,Lig2,-3.45,0.12\n")
        f.write("Lig2,Lig3,1.85,0.10\n")
        f.write("Lig3,Lig4,-0.75,0.14\n")
        f.write("Lig4,Lig1,2.28,0.11\n")
        f.write("Lig1,Lig3,-1.55,0.13\n")
    print(f"Generated sample network CSV at: {os.path.abspath(net_file)}")


if __name__ == "__main__":
    generate_sample_fep_data()
