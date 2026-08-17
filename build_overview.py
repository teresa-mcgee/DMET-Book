#!/usr/bin/env python3
"""Generator for DMET-results/overview.qmd -- a whole-panel summary page.

Companion to build_chapters.py (same ROOT/PROTEINS conventions), but reads
project-wide result files instead of per-protein figures/, since this page
summarizes across all 27 proteins at once:

  - out/SignificantResults.txt         (DMETPAPER root) -- protein pQTL hits,
                                         permutation/GEV-adjusted. Authoritative
                                         for "is this protein's peak genuinely
                                         genome-wide significant".
  - output/pQTL_pandq.txt              (DMETPAPER root) -- top-locus p-value +
                                         cross-protein q-value, ALL 27 proteins.
  - out/SignificantResults_90.txt      (DMETPAPER root) -- RNA-side p/q + combined
                                         q, only where the pipeline actually ran
                                         to completion (currently just Ent1 --
                                         see dev-notes/DEV_NOTES.md for the known
                                         i/j indexing bug in
                                         R/dev/final_figure_generation.R blocking
                                         wider coverage).
  - out/mediator_BF_summary.txt        (DMETPAPER root) -- per-gene log10 Bayes
                                         Factors from saved bmediatR RDS objects,
                                         built by R/build_overview_mediator_table.R.
                                         Only covers proteins with a saved RDS
                                         (see that script's header for which).
  - out/TIMBR_allelic_series_summary.txt (DMETPAPER root) -- top allelic series +
                                         top number-of-alleles per protein with a
                                         saved TIMBR posterior object. OPTIONAL:
                                         this file may not exist yet (TIMBR rerun
                                         for Cyp2c50/Ces2 in progress) -- the
                                         script degrades gracefully if so.

Run this after re-running any of the above upstream data extractions, same as
build_chapters.py after a figure regeneration.

Usage: python3 build_overview.py
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DMETPAPER = ROOT.parent

PROTEINS = [
    l.strip().strip('"')
    for l in (ROOT / "biomarkernames_pyth.txt").read_text().splitlines()
    if l.strip()
]


def read_table(path, sep=None, has_quotes=False):
    """Minimal whitespace/tab table reader -> list of dict rows. No pandas dependency.

    Header is always split on generic whitespace (some files, e.g.
    SignificantResults_90.txt, have a space-separated header but tab-separated
    data rows -- an upstream inconsistency, not a typo here). Data rows split on
    `sep`. Also handles R's default write.table() leading unnamed row-index
    column (one more field than header) by dropping it."""
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    if not lines:
        return []
    header = lines[0].split()
    if has_quotes:
        header = [h.strip('"') for h in header]
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split(sep) if sep else line.split()
        if has_quotes:
            parts = [p.strip('"') for p in parts]
        if len(parts) == len(header) + 1:
            parts = parts[1:]  # drop R's unnamed row-index column
        if len(parts) != len(header):
            continue
        rows.append(dict(zip(header, parts)))
    return rows


def fmt_p(x, sig=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "—"
    if v != v:  # NaN
        return "—"
    if v == 0:
        return "0"
    return f"{v:.{sig}g}"


def load_protein_significance():
    """protein -> dict(sig_chrs=[...], best_row=dict) from out/SignificantResults.txt.
    A protein can have >1 significant chromosome (e.g. Ces2); best_row picks the
    lowest raw P.Val as the representative peak, matching build_chapters.py's
    scan_active_index() choice for the carousel, for consistency across the book."""
    rows = read_table(DMETPAPER / "out" / "SignificantResults.txt", sep=None)
    by_protein = {}
    for r in rows:
        p = r["Protein"]
        by_protein.setdefault(p, []).append(r)
    result = {}
    for p, rs in by_protein.items():
        rs_sorted = sorted(rs, key=lambda r: float(r["P.Val"]))
        result[p] = {"rows": rs_sorted, "best": rs_sorted[0]}
    return result


def load_pandq():
    """protein -> dict(pvalue, qvalue) from output/pQTL_pandq.txt. Covers all 27."""
    rows = read_table(DMETPAPER / "output" / "pQTL_pandq.txt", sep=" ", has_quotes=True)
    return {r["PROTEIN"]: r for r in rows}


def load_rna_sig90():
    """protein -> row from out/SignificantResults_90.txt. Sparse (see module docstring)."""
    rows = read_table(DMETPAPER / "out" / "SignificantResults_90.txt", sep="\t")
    return {r["Protein"]: r for r in rows}


def load_mediators():
    """List of dict rows from out/mediator_BF_summary.txt, or [] if not yet built."""
    return read_table(DMETPAPER / "out" / "mediator_BF_summary.txt", sep="\t")


def load_timbr():
    """protein -> row from out/TIMBR_allelic_series_summary.txt, or {} if pending."""
    rows = read_table(DMETPAPER / "out" / "TIMBR_allelic_series_summary.txt", sep="\t")
    return {r["Protein"]: r for r in rows}


def build_summary_table(sig, pandq, rna90, timbr):
    header = ("| Protein | pQTL Sig. | P-value | Adj. P-value | Q-value | "
              "RNA Sig. | RNA P-value | RNA Q-value | TIMBR Top Allelic Series | TIMBR # Alleles |")
    sep = "|---|---|---|---|---|---|---|---|---|---|"
    rows = [header, sep]
    for protein in PROTEINS:
        s = sig.get(protein)
        pq = pandq.get(protein, {})
        r90 = rna90.get(protein)
        tb = timbr.get(protein)

        pqtl_sig = "✅" if s else "❌"
        pval = fmt_p(pq.get("Pvalue"))
        adj_pval = fmt_p(s["best"]["Adjusted_P.val"]) if s else "—"
        qval = fmt_p(pq.get("Qvalue"))

        if r90:
            rna_sig = "✅"
            rna_pval = fmt_p(r90.get("P.Val"))
            rna_qval = fmt_p(r90.get("qvalue_RNA"))
        else:
            rna_sig = "—"
            rna_pval = "—"
            rna_qval = "—"

        if tb:
            top_series = tb.get("TopAllelicSeries", "—")
            top_k = tb.get("TopNumAlleles", "—")
        else:
            top_series = "_pending_"
            top_k = "_pending_"

        rows.append(
            f"| {protein} | {pqtl_sig} | {pval} | {adj_pval} | {qval} "
            f"| {rna_sig} | {rna_pval} | {rna_qval} | {top_series} | {top_k} |"
        )
    return "\n".join(rows)


def build_mediator_table(mediators):
    passing = [m for m in mediators
               if m.get("PassesThreshold") == "TRUE" and m.get("Gene.Symbol") != m.get("Protein")]
    if not passing:
        return "_No mediators currently exceed their threshold (excluding self-pairs)._"
    passing.sort(key=lambda m: (m["Protein"], m["MediatorType"], -float(m["LogBF"])))
    header = "| Protein | Type | Mediator | Log10 BF | Threshold | Threshold Source |"
    sep = "|---|---|---|---|---|---|"
    rows = [header, sep]
    for m in passing:
        kind = "RNA transcript" if m["MediatorType"] == "rna" else "Protein"
        rows.append(
            f"| {m['Protein']} | {kind} | {m['Gene.Symbol']} | {float(m['LogBF']):.2f} "
            f"| {float(m['ThresholdLogBF']):.2f} | {m['ThresholdSource']} |"
        )
    return "\n".join(rows)


def main():
    sig = load_protein_significance()
    pandq = load_pandq()
    rna90 = load_rna_sig90()
    mediators = load_mediators()
    timbr = load_timbr()

    L = []
    L += ['# Overview {.unnumbered}', '']
    L += ['This page summarizes pQTL/eQTL significance, mediation, and TIMBR '
          'allelic-series results across all 27 proteins in one place. Each '
          'section notes its data source and current coverage -- some columns '
          'are populated for every protein, others only for the subset that has '
          'been carried through that analysis stage so far (see '
          '[Analysis Status](status.qmd) for the per-stage completion table).', '']

    L += ['## Significance Summary', '']
    L += ['pQTL significance (✅/❌) is the permutation/GEV-corrected genome-wide '
          'call from `out/SignificantResults.txt` -- the authoritative source, not '
          'the looser uncorrected per-protein files. P-value and Q-value come from '
          '`output/pQTL_pandq.txt` and cover all 27 proteins; Adj. P-value is only '
          'available for the proteins with a genuine hit. RNA columns are marked '
          '"—" (not "❌") where the RNA-side analysis has not been run through to '
          'a comparable significance call yet -- **not** evidence of a null result. '
          'TIMBR columns show "_pending_" for proteins whose posterior model object '
          'has not been (re-)computed yet.', '']
    L += [build_summary_table(sig, pandq, rna90, timbr), '']

    L += ['## Significant Mediators', '']
    L += ['Genes/proteins whose mediation Bayes Factor exceeds the genome-wide '
          'permutation threshold for that protein (`permutation` source), or '
          'exceeds BF > 1 as a default cutoff where no permutation threshold has '
          'been computed yet (`default_BF1` source -- treat these as provisional). '
          'Self-pairs (a protein mediating its own QTL) are excluded as '
          'uninformative. Coverage is currently limited to proteins with a saved '
          'bmediatR result object -- see `dev-notes/DEV_NOTES.md`.', '']
    L += [build_mediator_table(mediators), '']

    (ROOT / "overview.qmd").write_text("\n".join(L) + "\n")
    print("wrote overview.qmd")
    print(f"  protein significance: {len(sig)}/27 proteins with a genuine hit")
    print(f"  p/q-value coverage: {len(pandq)}/27 proteins")
    print(f"  RNA significance coverage: {len(rna90)}/27 proteins")
    print(f"  mediator rows loaded: {len(mediators)} (passing threshold: "
          f"{sum(1 for m in mediators if m.get('PassesThreshold') == 'TRUE')})")
    print(f"  TIMBR coverage: {len(timbr)}/27 proteins" if timbr else "  TIMBR coverage: pending (summary file not found)")


if __name__ == "__main__":
    main()
