"""
Command Line Interface (CLI) for FEPCert.
"""

import sys
import os
import argparse
import numpy as np

from fepcert import __version__
from fepcert.parsers.gromacs_dhdl import parse_gromacs_dhdl_directory, parse_gromacs_dhdl
from fepcert.parsers.generic_fep import parse_generic_fep_csv
from fepcert.parsers.network_csv import parse_perturbation_network_csv
from fepcert.core.scoring import assess_fep_quality
from fepcert.reporters.plot_generator import generate_fepcert_figures
from fepcert.reporters.manuscript_prep import generate_fepcert_manuscript_assets
from fepcert.reporters.html_report import generate_fepcert_html_report


def print_banner():
    banner = rf"""
  ______ ______ _____   _____          _   
 |  ____|  ____|  __ \ / ____|        | |  
 | |__  | |__  | |__) | |     ___ _ __| |_ 
 |  __| |  __| |  ___/| |    / _ \ '__| __|
 | |    | |____| |    | |___|  __/ |  | |_ 
 |_|    |______|_|     \_____\___|_|   \__| v{__version__}

 Alchemical Free Energy Quality-Control & Reproducibility Toolkit
 Monreal-Hernández et al., 2026
"""
    print(banner)


def run_demo(output_dir: str = "fepcert_demo_output"):
    """
    Executes a benchmark demonstration evaluating an alchemical transformation across 11 lambda states
    with phase space overlap matrix, TI integration, time convergence, and thermodynamic cycle closure.
    """
    print(f"\n[FEPCert] Running demonstration benchmark on Alchemical Transformation (LigA -> LigB)...")
    os.makedirs(output_dir, exist_ok=True)
    
    metadata = {
        "transformation": "Ligand A -> Ligand B (p38 MAP Kinase)",
        "engine": "GROMACS 2024.1",
        "temperature_k": 298.15
    }
    
    # 1. Generate 11 lambda windows (0.0 to 1.0) with realistic non-linear dH/dlambda gradients
    lambdas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    n_frames = 1000
    rng = np.random.default_rng(42)
    
    gradients = []
    for lam in lambdas:
        # Realistic S-curve gradient: -15 * lam^2 + 8 * lam - 2 + noise
        base_grad = -12.0 * (lam ** 1.5) + 6.0 * lam - 2.5
        noise = rng.normal(0.0, 1.2, size=n_frames)
        # Small stationary drift
        time_series = base_grad + noise
        gradients.append(time_series)

    # 2. Mock 4-ligand perturbation network (L1 -> L2 -> L3 -> L4 -> L1)
    network_edges = [
        {"ligand_a": "Lig1", "ligand_b": "Lig2", "delta_g": -3.45, "error": 0.12},
        {"ligand_a": "Lig2", "ligand_b": "Lig3", "delta_g": 1.85, "error": 0.10},
        {"ligand_a": "Lig3", "ligand_b": "Lig4", "delta_g": -0.75, "error": 0.14},
        {"ligand_a": "Lig4", "ligand_b": "Lig1", "delta_g": 2.28, "error": 0.11},  # Sum = -3.45 + 1.85 - 0.75 + 2.28 = -0.07 kcal/mol
        {"ligand_a": "Lig1", "ligand_b": "Lig3", "delta_g": -1.55, "error": 0.13}   # Triangle 1-2-3 closure: -3.45 + 1.85 - (-1.55) = -0.05
    ]

    print("  -> Performing Thermodynamic Integration, phase space overlap calculation, and hysteresis audit...")
    report = assess_fep_quality(
        metadata=metadata,
        lambda_values=lambdas,
        gradients_list=gradients,
        network_edges=network_edges,
        temperature_k=298.15,
        unit="kcal/mol"
    )

    print("  -> Generating publication-ready vector figures (overlap matrix heatmap, TI curve, and time-convergence)...")
    generate_fepcert_figures(report, output_dir)

    print("  -> Drafting manuscript Methods text snippet, summary LaTeX tables, and BibTeX citations...")
    assets = generate_fepcert_manuscript_assets(report, output_dir)

    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()

    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing interactive report to {html_p}...")
    generate_fepcert_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)

    print("\n" + "="*70)
    print(f" [RESULT] Overall Alchemical Certification Status: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    print(f" * Transformation   : {report.metadata['transformation']}")
    print(f" * Free Energy (TI) : Delta G = {report.free_energy.delta_g:.2f} \u00b1 {report.free_energy.delta_g_error:.2f} {report.free_energy.unit} | Status: {report.free_energy.status}")
    if report.overlap_result:
        print(f" * Lambda Overlap   : Min adjacent Pi = {report.overlap_result.min_adjacent_overlap*100:.1f}% (Mean: {report.overlap_result.mean_adjacent_overlap*100:.1f}%) | Status: {report.overlap_result.status}")
    if report.convergence_result:
        print(f" * Time Convergence : Hysteresis = {report.convergence_result.final_hysteresis:.2f} kcal/mol (W_diss = {report.convergence_result.dissipated_work:.2f}) | Status: {report.convergence_result.status}")
    if report.cycle_result:
        print(f" * Cycle Closure    : RMSE = {report.cycle_result.cycle_rmse:.2f} kcal/mol ({report.cycle_result.n_cycles} cycles evaluated) | Status: {report.cycle_result.status}")
    print("="*70)
    print(f"\nAll outputs successfully saved to: {os.path.abspath(output_dir)}/")
    print(f"Open {os.path.abspath(html_p)} in your browser to inspect the full report.\n")


def run_assess(args):
    """
    Evaluates user-provided GROMACS directory or CSV file.
    """
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    if args.dir:
        print(f"\n[FEPCert] Parsing GROMACS dhdl.xvg files from directory: {args.dir}...")
        data = parse_gromacs_dhdl_directory(args.dir, file_pattern=args.pattern or "dhdl*.xvg")
        lambdas = data["lambda_values"]
        grads = data["gradients_list"]
    elif args.input:
        print(f"\n[FEPCert] Parsing alchemical table from: {args.input}...")
        data = parse_generic_fep_csv(args.input)
        lambdas = data["lambda_values"]
        grads = data["gradients_list"]
    else:
        print("[Error] Please specify either an input directory with --dir or a file with --input.", file=sys.stderr)
        sys.exit(1)

    network_edges = None
    if args.network:
        print(f"  -> Parsing perturbation network from: {args.network}...")
        network_edges = parse_perturbation_network_csv(args.network)

    meta = {
        "transformation": args.name or "Alchemical Transformation",
        "engine": args.engine or "GROMACS",
        "temperature_k": float(args.temperature)
    }

    print("  -> Performing free energy estimation and convergence audit...")
    report = assess_fep_quality(
        metadata=meta,
        lambda_values=lambdas,
        gradients_list=grads,
        network_edges=network_edges,
        temperature_k=float(args.temperature),
        unit=args.unit
    )

    print("  -> Generating publication figures...")
    generate_fepcert_figures(report, output_dir)

    print("  -> Generating manuscript text, LaTeX summary table, and BibTeX citations...")
    assets = generate_fepcert_manuscript_assets(report, output_dir)

    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()

    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing HTML quality report to {html_p}...")
    generate_fepcert_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)

    print("\n" + "="*70)
    print(f" [RESULT] Overall Quality Certification: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    print(f" * Free Energy : Delta G = {report.free_energy.delta_g:.2f} \u00b1 {report.free_energy.delta_g_error:.2f} {report.free_energy.unit}")
    if report.overlap_result:
        print(f" * Overlap     : Min adjacent Pi = {report.overlap_result.min_adjacent_overlap*100:.1f}%")
    print("="*70)
    print(f"\nReport ready at: {os.path.abspath(html_p)}\n")


def run_cycle(args):
    """
    Evaluates thermodynamic cycle closure directly from a network CSV file.
    """
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n[FEPCert] Parsing perturbation network from: {args.input}...")
    edges = parse_perturbation_network_csv(args.input)
    
    from fepcert.core.cycles import evaluate_cycle_closure
    cycle_res = evaluate_cycle_closure(edges)
    
    print("\n" + "="*70)
    print(f" [RESULT] Thermodynamic Cycle Closure Status: {cycle_res.status}")
    print("="*70)
    print(f" * Nodes (Ligands) : {cycle_res.n_nodes}")
    print(f" * Edges (Runs)    : {cycle_res.n_edges}")
    print(f" * Closed Cycles   : {cycle_res.n_cycles}")
    print(f" * Cycle RMSE      : {cycle_res.cycle_rmse:.3f} kcal/mol")
    print(f" * Max Cycle Error : {cycle_res.max_cycle_error:.3f} kcal/mol")
    print("="*70)
    print(f"Diagnostic: {cycle_res.diagnostic_message}\n")


def print_citation():
    bib = """@software{monreal2026fepcert,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{FEPCert: An Open-Source Toolkit for Quality-Control, Phase Space Overlap, and Thermodynamic Cycle Closure Certification of Alchemical Free Energy Simulations}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/fepcert}
}"""
    print("\nIf you use FEPCert in your publications, please cite:\n")
    print("APA Style:")
    print("Monreal-Hernández, A. (2026). FEPCert: An Open-Source Toolkit for Quality-Control, Phase Space Overlap, and Thermodynamic Cycle Closure Certification of Alchemical Free Energy Simulations (v1.0.0). Zenodo. https://github.com/sircalch/fepcert\n")
    print("BibTeX:")
    print(bib)
    print()


def main():
    parser = argparse.ArgumentParser(
        prog="fepcert",
        description="FEPCert: Quality-Control, Phase Space Overlap, and Thermodynamic Cycle Closure for Alchemical Free Energy Simulations."
    )
    parser.add_argument("-v", "--version", action="version", version=f"fepcert {__version__}")

    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # Assess command
    assess_parser = subparsers.add_parser("assess", help="Assess alchemical free energy calculations (TI / BAR / Overlap)")
    assess_parser.add_argument("-i", "--input", default=None, help="Path to alchemical CSV/TSV table")
    assess_parser.add_argument("-d", "--dir", default=None, help="Directory containing GROMACS dhdl*.xvg files")
    assess_parser.add_argument("--pattern", default="dhdl*.xvg", help="File pattern for GROMACS files (default: dhdl*.xvg)")
    assess_parser.add_argument("--network", default=None, help="Path to perturbation network CSV for cycle closure")
    assess_parser.add_argument("-o", "--output", default="fepcert_output", help="Directory for output report and assets (default: fepcert_output)")
    assess_parser.add_argument("--name", default=None, help="Transformation description (e.g. 'Lig1 -> Lig2')")
    assess_parser.add_argument("--engine", default="GROMACS", help="Simulation engine (GROMACS, AMBER, NAMD, OpenMM)")
    assess_parser.add_argument("--temperature", default=298.15, help="Simulation temperature in Kelvin (default: 298.15)")
    assess_parser.add_argument("--unit", default="kcal/mol", help="Energy unit: kcal/mol or kJ/mol (default: kcal/mol)")

    # Cycle command
    cycle_parser = subparsers.add_parser("cycle", help="Audit thermodynamic cycle closure in chemical perturbation network")
    cycle_parser.add_argument("-i", "--input", required=True, help="Path to network CSV file")
    cycle_parser.add_argument("-o", "--output", default="fepcert_cycle_output", help="Output directory (default: fepcert_cycle_output)")

    # Demo command
    demo_parser = subparsers.add_parser("demo", help="Run FEPCert on an alchemical transformation benchmark (11 lambda windows + 4-ligand cycle)")
    demo_parser.add_argument("-o", "--output", default="fepcert_demo_output", help="Output directory (default: fepcert_demo_output)")

    # Cite command
    subparsers.add_parser("cite", help="Display BibTeX and APA citation details")

    if len(sys.argv) == 1:
        print_banner()
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    if args.command == "assess":
        print_banner()
        run_assess(args)
    elif args.command == "cycle":
        print_banner()
        run_cycle(args)
    elif args.command == "demo":
        print_banner()
        run_demo(args.output)
    elif args.command == "cite":
        print_banner()
        print_citation()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

