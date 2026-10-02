#!/usr/bin/env python3
"""Copy v2 significant-locus analyses and build a locus-specific book page."""
import csv
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
TAG = "rerun_2026-09-24_v2"
SOURCE = PROJECT / "out" / TAG
SOURCE_FIG = SOURCE / "figs"
FIG_ROOT = ROOT / "figures" / TAG / "credible_interval_loci"


def read_significant_rows():
    with (SOURCE / "SignificantResults.txt").open() as f:
        rows = list(csv.DictReader(f, delimiter=" "))
    rows = [r for r in rows if r["Protein"].lower() != "ces2"]
    return rows


def copy_asset(source, folder, dest_name):
    if not source.exists():
        return None
    folder.mkdir(parents=True, exist_ok=True)
    dest = folder / dest_name
    shutil.copy2(source, dest)
    return dest.relative_to(ROOT).as_posix()


def link(path, label):
    return f"[{label}]({path})" if path else "Unavailable in this run"


def main():
    rows = read_significant_rows()
    FIG_ROOT.mkdir(parents=True, exist_ok=True)
    page = [
        "# Significant pQTL credible-interval follow-up {.unnumbered}",
        "",
        f"This page follows the significant loci in `{TAG}/SignificantResults.txt`. "
        "It uses the saved 90% bootstrap credible intervals and matching v2 figures. "
        "CES2 is omitted because of the protein QC issue. Each locus is linked to its "
        "genome scan, CI/gene plot, mediation, TIMBR CI sweep, and Merge result where saved.",
        "",
        "All coordinates are Mb. The listed interval is the 90% bootstrap CI. "
        "Mediation, TIMBR, and Merge links refer to the same protein and lead locus.",
        "",
        "| Protein | Chr | Lead locus | Mb | 90% CI (Mb) | Scan p | Adjusted p |",
        "|---|---:|---|---:|---:|---:|---:|",
    ]
    for r in rows:
        protein, chrom, locus = r["Protein"], r["CH"], r["Locus"]
        protein_label = protein.upper()
        lo, hi = float(r["CI.Start"]), float(r["CI.End"])
        page.append(
            f"| {protein_label} | {chrom} | {locus} | {float(r['Mb']):.2f} | "
            f"{lo:.2f}–{hi:.2f} | {float(r['P.Val']):.2e} | "
            f"{float(r['Adjusted_P.val']):.2e} |"
        )

    missing = []
    for r in rows:
        protein, chrom, locus = r["Protein"], str(r["CH"]), r["Locus"]
        protein_label = protein.upper()
        lo, hi = float(r["CI.Start"]), float(r["CI.End"])
        locus_dir = FIG_ROOT / f"{protein}_chr{chrom}_{locus}"
        src_full = SOURCE_FIG / protein / "full"
        src_timbr = SOURCE_FIG / protein / "TIMBR_ci"
        src_data = SOURCE / protein / "full"

        ci_plot = copy_asset(src_full / f"ci_genes_chr{chrom}.png", locus_dir,
                             "credible_interval_gene_track.png")
        scan_pdf = copy_asset(src_full / "scans_and_adjusted0.95.pdf", locus_dir,
                              "genome_scan_v2.pdf")
        med_overlay = copy_asset(src_full / f"mediation_overlay_chr{chrom}_withgenes.png",
                                 locus_dir, "mediation_overlay_withgenes.png")
        med_post_nogenes = copy_asset(src_full / f"mediation_overlay_chr{chrom}_nogenes.png",
                                      locus_dir, "mediation_posterior_nogenes.png")
        med_post_withgenes = med_overlay
        med_bf_nogenes = copy_asset(src_full / f"mediation_overlay_bf_chr{chrom}_nogenes.png",
                                    locus_dir, "mediation_bf_nogenes.png")
        med_bf_withgenes = copy_asset(src_full / f"mediation_overlay_bf_chr{chrom}_withgenes.png",
                                      locus_dir, "mediation_bf_withgenes.png")
        med_bars = copy_asset(src_full / f"targeted_mediation_chr{chrom}_posterior_bars.pdf",
                              locus_dir, "mediation_posterior_bars.pdf")
        timbr_dot = copy_asset(src_timbr / f"timbr_ci_dot_chr{chrom}.jpg", locus_dir,
                               "timbr_ci_sweep.jpg")
        timbr_bf = copy_asset(src_timbr / f"timbr_ci_bf_dot_chr{chrom}.jpg", locus_dir,
                              "timbr_ci_bf.jpg")
        merge = copy_asset(src_full / f"mergeplot_{locus}.png", locus_dir,
                           "merge_at_locus.png")
        ci_text = copy_asset(src_data / f"{locus}_ci.txt", locus_dir,
                             "credible_interval_scan_data.txt")
        haplotypes = []
        for f in sorted(src_timbr.glob(f"haplotype_top*_chr{chrom}.jpg")):
            haplotypes.append(copy_asset(f, locus_dir, f.name))

        page.extend([
            "",
            f"## {protein_label} chr{chrom}: {locus}",
            "",
            f"Lead position: {float(r['Mb']):.2f} Mb; 90% CI: {lo:.2f}–{hi:.2f} Mb; "
            f"scan p = {float(r['P.Val']):.2e}; adjusted p = {float(r['Adjusted_P.val']):.2e}.",
            "",
        ])
        if ci_plot:
            page.extend([f"![90% bootstrap credible interval and nearby genes for {protein_label}, chr{chrom}.]({ci_plot})", ""])
        else:
            page.extend(["**90% CI plot:** unavailable in the saved v2 figures.", ""])
            missing.append(f"{protein} chr{chrom}: CI plot")
        links = [
            f"Genome scan: {link(scan_pdf, 'v2 scan and adjusted scan PDF')}",
            "Mediation overlays: " + "; ".join([
                link(med_post_nogenes, "posterior, no gene track"),
                link(med_post_withgenes, "posterior, with gene track"),
                link(med_bf_nogenes, "Bayes factor, no gene track"),
                link(med_bf_withgenes, "Bayes factor, with gene track"),
                link(med_bars, "posterior bars"),
            ]),
            f"TIMBR: {link(timbr_dot, 'posterior-probability plot')}; {link(timbr_bf, 'Bayes-factor plot')}",
            "Top TIMBR haplotype plots: " + (", ".join(f"[{Path(x).name}]({x})" for x in haplotypes) if haplotypes else "not saved"),
            f"Merge: {link(merge, 'locus-specific merge plot')}",
            f"CI scan data: {link(ci_text, 'download')}",
        ]
        page.extend(["- " + x for x in links])
        page.append("")
        for caption, asset in (
            ("Mediation posterior probabilities without gene track", med_post_nogenes),
            ("Mediation posterior probabilities with gene track", med_post_withgenes),
            ("Mediation Bayes factors without gene track", med_bf_nogenes),
            ("Mediation Bayes factors with gene track", med_bf_withgenes),
            ("TIMBR posterior probabilities", timbr_dot),
            ("TIMBR Bayes factors", timbr_bf),
            ("Merge profile for this locus and interval", merge),
        ):
            if asset:
                page.extend([f"![{caption}: {protein_label} chr{chrom}.]({asset})", ""])
        for label, path in (("mediation posterior", med_post_withgenes),
                            ("mediation Bayes factor", med_bf_withgenes),
                            ("TIMBR posterior", timbr_dot),
                            ("TIMBR Bayes factor", timbr_bf), ("Merge", merge)):
            if path is None:
                missing.append(f"{protein} chr{chrom}: {label}")

    page.extend([
        "## Coverage notes",
        "",
        f"The v2 significant-results file contains {len(rows)} significant non-CES2 loci. "
        "Each is shown above with its 90% interval. Merge figures are matched to the lead locus "
        "and interval; plots from other regions are not substituted.",
        "",
        "Mediation overlays are the saved v2 targeted results and should be interpreted as "
        "locus follow-up, not a genome-wide mediation scan. TIMBR plots summarize the "
        "credible-interval sweep and the linked haplotype plots show the top saved partitions.",
    ])
    (ROOT / "credible_intervals.qmd").write_text("\n".join(page) + "\n")
    print(f"Built credible_intervals.qmd for {len(rows)} significant non-CES2 loci")
    if missing:
        print("Unavailable locus-specific items:")
        for x in missing:
            print(" -", x)


if __name__ == "__main__":
    main()
