# DMET-results — Developer Notes

Personal reference only. Not part of the published book (not in `_quarto.yml`),
not overwritten by `build_chapters.py`. Update by hand as issues are found/fixed.
The reader-facing summary of most of this lives in `status.qmd`, but that file is
regenerated from scratch by `build_chapters.py` — anything hand-written there gets
clobbered on the next run, hence this file.

## Open bugs

### Per-protein heritability figure is broken (all 27 proteins)
Every chapter's heritability image is the *same* placeholder page pulled from
`individual_h2.pdf`. `individual_h2.pdf` does have real per-protein estimates
(27 pages, one per protein), but pages carry no protein label/title — only a
locus ID and h2 value — so there's no reliable way to recover the page→protein
mapping from the PDF alone. Would need to re-derive per-protein h2 from the
original scan/model objects instead of parsing the PDF.
Generator-side flag: `build_chapters.py` lines ~201, ~251-253.

### RNA-mediation image filename mismatch (worked around, not fixed)
Figure generator outputs `genome_wide_rna_med_page_001.jpg`; older chapter
template hardcoded `genome_wide_mediators_page_001.jpg` (only matched Ent1).
`build_chapters.py` now checks both names, so this no longer breaks anything,
but the two filenames diverging is itself a smell — worth unifying at the
generator level if `R/dev/final_figure_generation.R` is touched again.

## Content gaps (not bugs — analysis just hasn't been run)

- **Merge/diplotype fine-mapping**: only run for Bsep and Cyp2c39, plus a
  partial Cyp2c50 run (diplotype filtering done, no significant SNPs surfaced).
  ~24 proteins have no merge output at all. Scripts live in `R/Merge/`
  (`Bsep_Merge.R`, `cyp2c39merge.R`, `cyp_Merge.R`, `sub_cyp_Merge.R`) and would
  need generalizing to run on the rest of the panel.
- **pQTL scans**: several proteins missing box-cox/RINT/full-dataset scan
  variants — Oct1.2, Ugt1a1, Ugt1a6 have essentially no scan output yet.
- **RNA mediation**: about half the panel has no mediation analysis run;
  a few more have the RDS result but no rendered genome-wide figure.
- **Missing bibliography entries**: Qasem et al. 2020 (QTAP protocol) and
  Oreper et al. 2017 (Inbred Strain Variant DB) are cited by author/year only
  in `intro.qmd` — full entries never added to `references.bib`. Flagged
  inline in `references.bib` and `intro.qmd` rather than fabricated.

## Housekeeping / risks

- **Outer `DMETpaper` repo has never been committed to git** (confirmed
  2026-08-16 — `git status` there still shows "No commits yet" with everything
  untracked). `DMET-results/` is a separate, already-committed/pushed repo
  nested inside it — don't confuse the two when checking history or pushing.
- **`DESCRIPTION`/`NAMESPACE` at the outer repo root are unfilled boilerplate**
  — never customized, not a real installable package. `devtools::load_all()`
  etc. won't work.
- **Legacy analysis scripts still hardcode paths into archived figure trees**
  (`R/pqtl_all.R`, `R/run_zeroInflatedScans.R`, `R/zero_inf_prot.R`, most of
  `R/dev/`, `R/TIMBR/dev/`, `R/mediation/dev|complete/`) — they write to
  `figs/`/`figures_legacy/`, which are superseded and gitignored. If any of
  these are resurrected, redirect their output to
  `DMET-results/figures/{protein}/...` first.

## Fixed (kept here for history — see status.qmd "Known issues fixed" too)

- Duplicate/stale figure trees (`figs/`, top-level `figures/`, a nested
  `figures/figures/` copy) merged into the single `figures/{protein}/` tree
  that `run_final_fig.sh` writes to.
- TIMBR sections were `include=FALSE`/commented out in every chapter despite
  real output existing for 26/27 proteins (all but Ces2) — now wired in.
