#!/usr/bin/env python3
"""Refresh the v2 overview, progress page, and reporting figure bundle."""
import csv
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
TAG = "rerun_2026-09-24_v2"
RNA_TAG = "rerun_2026-09-27_rna_v2"
SOURCE = PROJECT / "out" / TAG
FIG_SOURCE = SOURCE / "figs" / "all"
FIG_DIR = ROOT / "figures" / TAG
RNA_FIG_DIR = ROOT / "figures" / RNA_TAG
MEDIATION_SOURCE = PROJECT / "out" / RNA_TAG / "figs"
MEDIATION_PROTEINS = [
    x.strip() for x in (PROJECT / "slurm" / "rerun2026" / "rna_v2_mediation_proteins.txt").read_text().splitlines()
    if x.strip() and x.strip().lower() != "ces2"
]
H2_CACHE = SOURCE / "coding_transcript_h2_broad_sense.csv"
CORR_FIG_DIR = ROOT / "figures" / TAG


def load_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def val(row, key):
    x = row.get(key, "")
    return "" if x in (None, "NA", "NaN", "nan") else x


def number(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def fmt(x, digits=2):
    v = number(x)
    return "—" if v is None else f"{v:.{digits}f}"


def fmtp(x):
    v = number(x)
    return "—" if v is None else f"{v:.2e}"


def md(x):
    return str(x).replace("|", "\\|").replace("\n", " ") if x else "—"


def main():
    subprocess.run([sys.executable, str(ROOT / "build_credible_interval_appendix.py")], check=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for f in FIG_SOURCE.iterdir():
        if f.suffix.lower() in {".png", ".pdf"}:
            shutil.copy2(f, FIG_DIR / f.name)
    for pattern in ("figure3_pqtl_scan_values.csv", "figure3_pqtl_significant_peaks.csv",
                    "figure3_pqtl_thresholds_90_percent.csv", "figure3_chr*_overlay_alpha.csv"):
        for f in FIG_SOURCE.glob(pattern):
            shutil.copy2(f, FIG_DIR / f.name)
    # Figure 2 and the companion clustering outputs are generated under the
    # established paper directory; copy their refreshed display versions into
    # the v2 reporting bundle after the run-level files above.
    paper_figures = ROOT / "figures" / "paper"
    for pattern in ("figure2_heritability_profiles.*", "protein_cluster_gapstat.*",
                    "protein_cluster_heatmap.*", "protein_pca_biplot*.*"):
        for f in paper_figures.glob(pattern):
            if f.suffix.lower() in {".png", ".pdf"}:
                shutil.copy2(f, FIG_DIR / f.name)
    for name in ("figure2_broad_sense_heritability_estimates.csv",
                 "protein_cluster_kw_drivers.csv", "protein_pca_loadings.csv",
                 "protein_pca_loadings_boxcox.csv", "protein_pca_loadings_zscore.csv"):
        f = paper_figures / name
        if f.exists():
            shutil.copy2(f, FIG_DIR / name)
    # Use the complete 2026-09-27 RNA eQTL map and publish the current targeted
    # mediation and TIMBR plot variants, omitting CES2 per the project QC
    # instruction. These plot-only outputs live under SOURCE/figs, while the
    # saved mediation analyses live under MEDIATION_SOURCE.
    RNA_FIG_DIR.mkdir(parents=True, exist_ok=True)
    mediation_dir = FIG_DIR / "mediation_rna_v2"
    mediation_dir.mkdir(exist_ok=True)
    timbr_dir = FIG_DIR / "timbr_ci_rna_v2"
    timbr_dir.mkdir(exist_ok=True)
    for protein in MEDIATION_PROTEINS:
        src = SOURCE / "figs" / protein / "full"
        bars_src = MEDIATION_SOURCE / protein / "full"
        dst = mediation_dir / protein
        if src.exists():
            dst.mkdir(exist_ok=True)
            for f in src.glob("mediation_overlay_*.png"):
                if f.stem.endswith("_example") or "_legacy_" in f.stem:
                    continue
                shutil.copy2(f, dst / f.name)
        if bars_src.exists():
            dst.mkdir(exist_ok=True)
            for f in bars_src.glob("targeted_mediation_*_posterior_bars.pdf"):
                shutil.copy2(f, dst / f.name)
        timbr_src = SOURCE / "figs" / protein / "TIMBR_ci"
        timbr_dst = timbr_dir / protein
        if timbr_src.exists():
            timbr_dst.mkdir(parents=True, exist_ok=True)
            for f in timbr_src.glob("timbr_ci_*.jpg"):
                shutil.copy2(f, timbr_dst / f.name)
    mediation_rows = []
    for protein in MEDIATION_PROTEINS:
        dst = mediation_dir / protein
        if not dst.exists():
            continue
        links = []
        for f in sorted(dst.glob("mediation_overlay_*.png")):
            if f.stem.endswith("_example") or "_legacy_" in f.stem:
                continue
            stem = f.stem.removeprefix("mediation_overlay_")
            stem = stem.replace("_nogenes", " no genes").replace("_withgenes", " with genes")
            stem = stem.replace("_bf_", " BF ").replace("_chr", " chr")
            links.append(f"[mediation {stem}](figures/{TAG}/mediation_rna_v2/{protein}/{f.name})")
        links.extend(
            f"[posterior bars chr{f.stem.split('_chr')[-1].split('_')[0]}](figures/{TAG}/mediation_rna_v2/{protein}/{f.name})"
            for f in sorted(dst.glob("targeted_mediation_*_posterior_bars.pdf"))
        )
        timbr_dst = timbr_dir / protein
        for f in sorted(timbr_dst.glob("timbr_ci_*.jpg")) if timbr_dst.exists() else []:
            variant = "Bayes factor" if "_bf_" in f.name else "posterior probability"
            chrom = f.stem.split("_chr")[-1]
            links.append(f"[TIMBR {variant} chr{chrom}](figures/{TAG}/timbr_ci_rna_v2/{protein}/{f.name})")
        mediation_rows.append(f"| {protein.upper()} | " + ", ".join(links) + " |")
    shutil.copy2(SOURCE / "coding_transcript_h2_broad_sense.csv",
                FIG_DIR / "coding_transcript_h2_broad_sense.csv")
    for ext in ("png", "pdf"):
        p = CORR_FIG_DIR / f"protein_phenotype_correlation_heatmap.{ext}"
        if not p.exists():
            raise FileNotFoundError(f"Create descriptive correlation plot first: {p}")

    rows = load_csv(SOURCE / "prot_summary_table_with_heritability.csv")
    rows = [r for r in rows if val(r, "Protein").upper() != "CES2"]
    # Preserve source-folder keys in their original case; publish the protein
    # symbols in uppercase in the summary tables.
    for r in rows:
        r["Protein"] = val(r, "Protein").upper()
    # The v2 pQTL summary omits Oct1.2, although its Box-Cox protein H2
    # estimate is present in the current broad-sense H2 results table. Retain
    # it in the integrated report with QTL fields explicitly left unavailable.
    if not any(val(r, "Protein").upper() == "OCT1.2" for r in rows):
        h2_path = ROOT / "figures" / "paper" / "figure2_broad_sense_heritability_estimates.csv"
        oct_h2 = next(r for r in load_csv(h2_path) if r["Protein"].upper() == "OCT1.2")
        oct_row = dict.fromkeys(rows[0].keys(), "")
        oct_row.update({
            "Protein": "OCT1.2", "Set": "full", "H2": oct_h2["H2"],
            "H2_CI_low": oct_h2["CI_low"], "H2_CI_high": oct_h2["CI_high"],
            "H2_Pvalue": oct_h2["p_value"],
            "H2_Pvalue_Bonferroni": oct_h2["p_value_bonferroni"],
            "Transcript_mapping": "No original-array gene-specific transcript match",
        })
        rows.append(oct_row)
    panel_order = [x.strip().upper() for x in (ROOT / "biomarkernames_pyth.txt").read_text().splitlines() if x.strip()]
    order = {protein: i for i, protein in enumerate(panel_order)}
    rows.sort(key=lambda r: order.get(val(r, "Protein").upper(), len(order)))
    # Publish a QC-filtered copy for the book; the original v2 source stays intact.
    report_csv = FIG_DIR / "v2_integrated_summary_no_CES2.csv"
    with report_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Individual v2 pQTL scan PDFs supplement the cross-panel plots.
    scan_dir = FIG_DIR / "pqtl_per_protein"
    scan_dir.mkdir(exist_ok=True)
    scan_files = {}
    for protein_fig_dir in (SOURCE / "figs").iterdir():
        scan = protein_fig_dir / "full" / "scans_and_adjusted0.95.pdf"
        if protein_fig_dir.is_dir() and scan.exists() and protein_fig_dir.name.upper() != "CES2":
            target = scan_dir / f"{protein_fig_dir.name}.pdf"
            shutil.copy2(scan, target)
            scan_files[protein_fig_dir.name.upper()] = target.name
    subprocess.run([sys.executable, str(ROOT / "update_protein_chapter_scans.py")], check=True)

    table = [
        "| Protein | pQTL lead (chr:Mb) | pQTL P | Adjusted P | Scan plots | Protein H² (95% CI) | Protein H² Bonferroni P | Coding transcript | Transcript H² (95% CI) | Transcript H² BH P |",
        "|---|---:|---:|---:|---|---:|---:|---|---:|---:|",
    ]
    for r in rows:
        locus = f"{val(r, 'CH')}:{fmt(val(r, 'Mb'), 2)}" if val(r, "CH") else "—"
        h2 = f"{fmt(val(r, 'H2'))} ({fmt(val(r, 'H2_CI_low'))}, {fmt(val(r, 'H2_CI_high'))})" if val(r, "H2") else "—"
        th2 = f"{fmt(val(r, 'Transcript_H2'))} ({fmt(val(r, 'Transcript_H2_CI_low'))}, {fmt(val(r, 'Transcript_H2_CI_high'))})" if val(r, "Transcript_H2") else "—"
        scan_link = f"[PDF](figures/{TAG}/pqtl_per_protein/{scan_files[val(r, 'Protein').upper()]})" if val(r, "Protein").upper() in scan_files else "—"
        table.append("| " + " | ".join([
            md(val(r, "Protein")), locus, fmtp(val(r, "P.Val")), fmtp(val(r, "Adjusted_P.val")),
            scan_link, h2, fmtp(val(r, "H2_Pvalue_Bonferroni")), md(val(r, "Transcript")), th2,
            fmtp(val(r, "Transcript_H2_Pvalue_BH")),
        ]) + " |")

    overview = f'''# V2 results overview {{.unnumbered}}

This is the reporting overview for the **{TAG}** rerun. It presents the current cross-protein results and links the figure set. Protein symbols are shown in uppercase. The reporting table and figures exclude **CES2** because of the previously identified protein QC issue. CES2's legacy source chapter is omitted from the rendered book; this page is the cross-panel summary.

## Cross-protein summary

The table combines pQTL lead-locus statistics and protein broad-sense H² with coding-transcript broad-sense H² where an original-array probe could be matched. Protein H² uses the Box–Cox transformed protein measurements supplied to the QTL analyses. Transcript H² uses RINT-transformed individual RNA array measurements, with CC line as a random effect and scan date adjusted as a fixed effect. Transcript rows without a gene-specific match are shown as unavailable.

Download the [QC-filtered CSV summary](figures/{TAG}/v2_integrated_summary_no_CES2.csv) or the [transcript H² estimates and intervals](figures/{TAG}/coding_transcript_h2_broad_sense.csv). OCT1.2 is retained with its protein H² estimate from the figure 2 results, but the v2 integrated table contains no OCT1.2 pQTL row.

{chr(10).join(table)}

## Heritability and strain structure

![Protein broad-sense H² and clustered Box–Cox protein profiles. Strains are ordered by the gap-statistic k=3 assignment.](figures/{TAG}/figure2_heritability_profiles.png)

Protein H² is estimated from replicated protein measurements. The figure displays the protein-level estimates and their uncertainty alongside the strain-level expression heatmap. The separately estimated transcript H² values are summarized in the table above; they are not the estimates in this forest panel.

### Coding-transcript heritability

![Broad-sense heritability estimates for the coding transcripts mapped to the protein panel.](figures/{TAG}/transcript_broad_sense_heritability_forest.png)

Transcript H² was estimated from RINT-transformed individual RNA measurements for 19 unique original-array probes, corresponding to 20 protein rows. The model includes a random CC-line effect and adjusts for scan date; intervals are 95% parametric-bootstrap intervals. The table gives the per-transcript BH-adjusted p-values.

### Protein phenotype correlation

![Descriptive across-strain Pearson correlations between protein abundance measurements.](figures/{TAG}/protein_phenotype_correlation_heatmap.png)

This descriptive Pearson matrix summarizes observed across-strain protein abundance and is distinct from the random-effect covariance estimates below.

### Genetic covariance and genetic correlation

![Protein genetic correlations, with broad-sense H² on the diagonal.](figures/{TAG}/protein_genetic_correlation_H2_diagonal.png)

![Coding-transcript genetic correlations, with broad-sense H² on the diagonal.](figures/{TAG}/transcript_genetic_correlation_H2_diagonal.png)

Download the [protein genetic covariance matrix](figures/{TAG}/protein_genetic_correlation_H2_diagonal_genetic_covariance.csv), [protein genetic correlation matrix with H² diagonal](figures/{TAG}/protein_genetic_correlation_H2_diagonal_genetic_correlation_H2_diagonal.csv), [coding-transcript covariance matrix](figures/{TAG}/transcript_genetic_correlation_H2_diagonal_genetic_covariance.csv), or [coding-transcript correlation matrix with H² diagonal](figures/{TAG}/transcript_genetic_correlation_H2_diagonal_genetic_correlation_H2_diagonal.csv). The displayed diagonal contains the existing broad-sense H² estimates; off-diagonal cells are estimated CC-line genetic correlations. Covariances use Box–Cox protein values and RINT transcript values. Protein residual cross-trait covariance is assumed zero because replicate assays are not paired across traits; transcript off-diagonal covariance is corrected for within-line residual covariance. Correlations were projected to a valid correlation matrix after sampling correction. Estimates are point estimates without confidence intervals.

## Genome-wide QTL overview

### pQTL

![Genome-wide pQTL scan signal across proteins.](figures/{TAG}/figure3_pqtl_genome_heatmap.png)

![Protein-specific pQTL signal overlaid on chromosome 3.](figures/{TAG}/figure3_pqtl_overlay_chr3.png)

![Protein-specific pQTL signal overlaid on chromosome 9.](figures/{TAG}/figure3_pqtl_overlay_chr9.png)

The v2 figure bundle also includes zoomed pQTL views for chromosomes 3, 16, 17, and 19: [chr3](figures/{TAG}/figure3_pqtl_zoom_chr3.png), [chr16](figures/{TAG}/figure3_pqtl_zoom_chr16.png), [chr17](figures/{TAG}/figure3_pqtl_zoom_chr17.png), and [chr19](figures/{TAG}/figure3_pqtl_zoom_chr19.png). The summary table links individual protein scan PDFs where available ({len(scan_files)} of {len(rows)} reportable proteins); current scan PDFs are missing for CYP2D9, CYP2D10, CYP2E1, MRP2, and OCT1.2.

### Significant-locus credible intervals

The [credible-interval appendix](credible_intervals.qmd) reports every significant non-CES2 v2 pQTL locus with its 90% bootstrap interval, genome scan, nearby-gene CI plot, and corresponding mediation, TIMBR, and Merge follow-up where available.

### eQTL

![eQTL genomic locations, with eQTL peak position on the x-axis and coding-gene position on the y-axis.](figures/{RNA_TAG}/figure3_eqtl_location_heatmap.png)

The merged corrected RNA-v2 eQTL table contains peak results for all 20,666 of 20,666 transcripts, with all 104 of 104 task outputs present. The location map plots the 19,959 transcripts with usable chromosome and genomic-position annotation on both axes; 707 eQTL results do not have complete positional annotation for this display. Color is capped at −log₁₀(p)=5, with dark outline circles marking peaks at or above this approximate significance threshold. The refreshed plot and its underlying coordinates are in the [RNA-v2 figure bundle](figures/{RNA_TAG}/).

## Protein variation and transcriptome context

### Protein PCA

![Protein PCA biplot.](figures/{TAG}/protein_pca_biplot.png)

Additional protein PCA versions are available for [Box–Cox values](figures/{TAG}/protein_pca_biplot_boxcox.png) and [z-scored values](figures/{TAG}/protein_pca_biplot_zscore.png). The [gap-statistic curve](figures/{TAG}/protein_cluster_gapstat.png) and [k=3 strain-cluster heatmap](figures/{TAG}/protein_cluster_heatmap.png) show the companion clustering analysis.

### Transcriptome

![RINT-normalized transcriptome hierarchical heatmap.](figures/{TAG}/transcriptome_rint_hierarchical_heatmap.png)

![Binned correlation structure across the RINT-normalized transcriptome.](figures/{TAG}/transcriptome_rint_all_transcript_binned_correlation.png)

The correlation overview bins transcripts by expression variance for display; it does not show a complete 20,666 × 20,666 transcript correlation matrix.

## RNA mediation coverage

The 2026-09-27 targeted mediation array was submitted for eight proteins. Seven are reportable after excluding CES2: {", ".join(p.upper() for p in MEDIATION_PROTEINS)}. The current posterior-probability and Bayes-factor overlays and TIMBR credible-region plots are copied into the book figure bundle below; this is targeted locus follow-up, not panel-wide mediation. Each overlay uses the Merge-style haplotype-association scan line and the original analysis points.

| Protein | Mediation and TIMBR figures |
|---|---|
{chr(10).join(mediation_rows)}

The report previously showed no mediation figures because its builder only copied cross-panel plots and pQTL scans from the 2026-09-24 folder; it did not collect the separate 2026-09-27 RNA mediation folder. The other 19 QC-retained proteins have no plots in this RNA-v2 folder because they were not included in the explicit SLURM mediation protein list. CES2 outputs are omitted from the report because of its protein QC issue.

## Remaining reporting gaps

- **Genetic covariance uncertainty:** current covariance and correlation matrices are point estimates; confidence intervals or bootstrap uncertainty are not yet estimated.
- **eQTL positional annotation:** 20,666/20,666 transcripts have peak results, while 19,959 are plotted because 707 lack complete position annotation.
- **Coding transcript matching:** six of the 26 QC-retained protein rows lack a gene-specific original-array transcript match, so their transcript H² is not estimated.
- **OCT1.2 pQTL:** its protein H² is available, but the v2 integrated summary and per-protein scan figure are absent.
- **Mediation and fine mapping:** mediation and TIMBR outputs are linked at all eight significant non-CES2 pQTL CIs. Merge plots are linked where available. These are locus-specific follow-up, not a uniform panel-wide report.
'''
    (ROOT / "overview.qmd").write_text(overview)

    status = f'''# Analysis progress {{.unnumbered}}

This page tracks the reporting status of the **{TAG}** results bundle. The cross-panel plots and integrated summary are in [V2 results overview](overview.qmd). CES2 is excluded from reporting analyses because of its protein QC issue.

| Analysis | Current status | Evidence and reporting note |
|---|---|---|
| Protein pQTL summary | Available for 25/26 QC-retained proteins | Lead locus, p-value, adjusted p-value, and thresholds are in the integrated summary; OCT1.2 has no v2 pQTL summary row. |
| Individual pQTL scan PDFs | Available for {len(scan_files)}/26 proteins | Linked from the integrated table where available; scan PDFs are missing for CYP2D9, CYP2D10, CYP2E1, MRP2, and OCT1.2. |
| Significant pQTL credible intervals | Available for 8 significant non-CES2 loci | The credible-interval page uses the 90% bootstrap CIs and v2 scan/CI figures. Mediation and TIMBR are linked for all 8; Merge figures are matched to each lead locus and interval. |
| Protein broad-sense H² | Available for all 26 QC-retained proteins | Box–Cox protein estimates with 95% bootstrap intervals are shown in the figure 2 panel; p-values include Bonferroni correction in the summary. |
| Protein strain profiles and clustering | Available | Box–Cox profiles, PCA, gap-statistic plot, and k=3 heatmap are linked from the overview. |
| Transcript broad-sense H² | Available for matched probes | Forest plot and summary table cover 19 unique probes mapped to 20 of 26 QC-retained protein rows. Per-probe estimates include 95% bootstrap intervals and adjusted p-values. |
| eQTL genome map | Peak results complete; location plot near-complete | The merged RNA-v2 table has 20,666/20,666 peak results and 104/104 task outputs; the location map plots 19,959 transcripts with complete positional annotation. |
| Transcriptome structure | Available as overview plots | RINT hierarchical heatmap and binned transcript-correlation display are included. |
| Genetic covariance / genetic correlations | Available as point estimates | Protein and coding-transcript CC-line covariance matrices and H²-diagonal correlation plots are included; uncertainty intervals have not been estimated. |
| Mediation | Partial, targeted | 7/26 reportable proteins have plots in the RNA-v2 targeted mediation folder; the submitted array explicitly listed eight proteins, including excluded CES2. |
| TIMBR at significant CIs | Available for all 8 significant non-CES2 loci | CI-sweep dot plots and saved top haplotype plots are linked from the credible-interval page. |
| Merge at significant CIs | Available for 7/8 significant non-CES2 loci | BSEP has no saved v2 Merge figure for its chr16 credible interval. |
| Merge/diplotype fine mapping | Partial | Available locus-specific results should be reported as selected-locus follow-up, not a complete panel-wide analysis. |

## Reporting interpretation

The completed plots support reporting protein-level heritability, strain variation, protein QTL signals, complete transcript peak results, near-complete positional eQTL mapping, and locus-specific CI follow-up. The current pQTL table is incomplete for OCT1.2. Genetic covariance estimates are available but currently lack uncertainty intervals; the BSEP significant-locus CI lacks a matching Merge analysis plot.
'''
    (ROOT / "status.qmd").write_text(status)
    print(f"Wrote overview.qmd and status.qmd; copied v2 figure bundle to {FIG_DIR}")


if __name__ == "__main__":
    main()
