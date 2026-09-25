"""
Quickstart API tutorial for FEPCert.
"""

import os
import sys

# Ensure current script dir is in sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fepcert import assess_fep_quality
from fepcert.parsers import parse_gromacs_dhdl_directory, parse_perturbation_network_csv
from fepcert.reporters import (
    generate_fepcert_figures,
    generate_fepcert_manuscript_assets,
    generate_fepcert_html_report
)
from generate_sample_dhdl_data import generate_sample_fep_data


def main():
    print("Running FEPCert Python API quickstart tutorial...")
    raw_data_dir = "sample_gromacs_dhdl"
    generate_sample_fep_data(raw_data_dir)
    
    output_dir = "quickstart_fepcert_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Parse simulation data and network
    fep_data = parse_gromacs_dhdl_directory(raw_data_dir, file_pattern="dhdl*.xvg")
    net_file = os.path.join(raw_data_dir, "sample_perturbation_network.csv")
    network_edges = parse_perturbation_network_csv(net_file)
    
    # 2. Assess free energy convergence
    report = assess_fep_quality(
        metadata={"transformation": "Lig1 -> Lig2", "engine": "SYNTHETIC DEMO DATA"},
        lambda_values=fep_data["lambda_values"],
        gradients_list=fep_data["gradients_list"],
        network_edges=network_edges,
        unit="kcal/mol"
    )
    
    print(f"\nOverall Certification: {report.overall_status}")
    print(f"Validation Score: {report.validation_score}")
    print(f"Free Energy (TI): Delta G = {report.free_energy.delta_g:.2f} \u00b1 {report.free_energy.delta_g_error:.2f} {report.free_energy.unit}")
    print(f"Min Adjacent Overlap: {report.overlap_result.min_adjacent_overlap*100:.1f}%")
    print(f"Cycle RMSE: {report.cycle_result.cycle_rmse:.2f} kcal/mol ({report.cycle_result.n_cycles} cycles)")
    
    # 3. Export all publication assets
    generate_fepcert_figures(report, output_dir)
    assets = generate_fepcert_manuscript_assets(report, output_dir)
    
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()
        
    html_p = os.path.join(output_dir, "report.html")
    generate_fepcert_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)
    
    print(f"\nCompleted! HTML report available at: {os.path.abspath(html_p)}")


if __name__ == "__main__":
    main()
