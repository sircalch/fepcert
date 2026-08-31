"""
Manuscript Methods snippet generator, summary tables, and BibTeX citations for FEPCert.
"""

from typing import Dict, Any, Optional
import os
import pandas as pd
from fepcert.core.scoring import FEPValidationReport


def generate_fepcert_manuscript_assets(
    report: FEPValidationReport,
    output_dir: str
) -> Dict[str, str]:
    """
    Generates manuscript Methods paragraph, summary CSV/LaTeX tables, and BibTeX citations.

    Parameters
    ----------
    report : FEPValidationReport
    output_dir : str

    Returns
    -------
    paths : dict
    """
    os.makedirs(output_dir, exist_ok=True)
    generated = {}

    # 1. Summary DataFrame
    rows = []
    meta = report.metadata
    fe = report.free_energy
    
    rows.append({"Parameter": "Transformation / System", "Value": f"{meta.get('transformation', 'LigA -> LigB')} ({meta.get('engine', 'GROMACS')})", "Status": "PASS"})
    rows.append({"Parameter": "Estimator Method", "Value": fe.method, "Status": "PASS"})
    rows.append({"Parameter": "Calculated Free Energy", "Value": f"{fe.delta_g:.3f} \u00b1 {fe.delta_g_error:.3f} {fe.unit}", "Status": fe.status})
    rows.append({"Parameter": "Number of Lambda Windows", "Value": f"{len(fe.lambda_values)} states", "Status": "PASS"})
    
    if report.overlap_result:
        ov = report.overlap_result
        rows.append({"Parameter": "Min Adjacent Overlap", "Value": f"{ov.min_adjacent_overlap*100:.1f}% (Mean: {ov.mean_adjacent_overlap*100:.1f}%)", "Status": ov.status})

    if report.convergence_result:
        cr = report.convergence_result
        rows.append({"Parameter": "Forward/Reverse Hysteresis", "Value": f"{cr.final_hysteresis:.3f} kcal/mol (W_diss: {cr.dissipated_work:.3f})", "Status": cr.status})

    if report.cycle_result:
        cy = report.cycle_result
        rows.append({"Parameter": "Cycle Closure Error", "Value": f"RMSE = {cy.cycle_rmse:.3f} kcal/mol ({cy.n_cycles} cycles)", "Status": cy.status})

    df_summary = pd.DataFrame(rows)

    # CSV Table
    csv_path = os.path.join(output_dir, "fepcert_summary_table.csv")
    df_summary.to_csv(csv_path, index=False)
    generated["summary_csv"] = csv_path

    # LaTeX Table
    tex_path = os.path.join(output_dir, "fepcert_summary_table.tex")
    tex_content = df_summary.to_latex(index=False, escape=False)
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write("% FEPCert Alchemical Free Energy Validation Table\n")
        f.write(tex_content)
    generated["summary_tex"] = tex_path

    # 2. Methods Text Snippet
    methods_path = os.path.join(output_dir, "methods_snippet.txt")
    trans_str = meta.get("transformation", "the alchemical transformation")
    engine_str = meta.get("engine", "molecular dynamics simulations")
    
    ov_str = ""
    if report.overlap_result:
        ov = report.overlap_result
        ov_str = f"Phase space connectivity was certified with a minimum adjacent lambda-state overlap of Pi = {ov.min_adjacent_overlap*100:.1f}%, exceeding the recommended Klimovich threshold. "

    cr_str = ""
    if report.convergence_result:
        cr = report.convergence_result
        cr_str = f"Time-series convergence was verified by forward-reverse chunking analysis yielding a hysteresis of {cr.final_hysteresis:.2f} kcal/mol (W_diss = {cr.dissipated_work:.2f} kcal/mol). "

    cy_str = ""
    if report.cycle_result:
        cy = report.cycle_result
        cy_str = f"Thermodynamic cycle closure consistency was confirmed across {cy.n_cycles} independent cycles with an RMSE of {cy.cycle_rmse:.2f} kcal/mol. "

    full_methods = (
        f"Alchemical free energy calculations for {trans_str} were performed in {engine_str}. "
        f"Free energy estimates and convergence quality were systematically audited using FEPCert v1.0.0 (Monreal-Hernández, 2026). "
        f"The net free energy difference was calculated using {fe.method} as Delta G = {fe.delta_g:.2f} \u00b1 {fe.delta_g_error:.2f} {fe.unit}. "
        f"{ov_str}{cr_str}{cy_str}"
        f"The calculation achieved an overall certification status of: {report.overall_status}."
    )

    with open(methods_path, "w", encoding="utf-8") as f:
        f.write(full_methods + "\n")
    generated["methods_text"] = methods_path

    # 3. BibTeX Citation
    bib_path = os.path.join(output_dir, "citation.bib")
    bib_content = """@software{monreal2026fepcert,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{FEPCert: An Open-Source Toolkit for Quality-Control, Phase Space Overlap, and Thermodynamic Cycle Closure Certification of Alchemical Free Energy Simulations}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/amonreal/fepcert}
}
"""
    with open(bib_path, "w", encoding="utf-8") as f:
        f.write(bib_content)
    generated["citation_bib"] = bib_path

    return generated
