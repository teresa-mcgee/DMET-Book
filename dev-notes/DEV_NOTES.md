# DMET-results — Developer Notes

Personal reference only. Not part of the published book (not in `_quarto.yml`),
not overwritten by `build_chapters.py`. Update by hand as issues are found/fixed.
The reader-facing summary of most of this lives in `status.qmd`, but that file is
regenerated from scratch by `build_chapters.py` — anything hand-written there gets
clobbered on the next run, hence this file.

## Open bugs

### Unreferenced duplicate `figures/{protein}h2.jpg` flat files
Alongside the real per-protein `figures/{protein}/h2.jpg`, there's a flat
`figures/{protein}h2.jpg` for all 27 proteins — these are all byte-identical
to each other (single placeholder page) and not referenced by any `.qmd`
chapter. Likely a leftover from the pre-dedup figure tree. Safe to delete;
not yet done because it's a batch delete across the tree — flagging instead
of doing it silently. `git ls-files figures/*h2.jpg` to list them.

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

- **Heritability figures were assumed broken but weren't — only ATPase and
  Abcb4 were swapped.** The 2026-08 cleanup wrongly concluded the whole
  `individual_h2.pdf` page→protein extraction never worked, based on the flat
  legacy duplicates (see "Open bugs" above) being byte-identical. The real
  `figures/{protein}/h2.jpg` files were correct for 25/27 proteins the whole
  time — only ATPase's and Abcb4's had been swapped (confirmed visually: each
  image's own title text named the other protein). Swapped the two files back
  2026-08-16; `build_chapters.py`'s `heritability` status now reports `"done"`
  instead of a hardcoded `"placeholder"`, and the per-chapter callout-warning
  about a broken pipeline was removed.
- **Genome-scan/TIMBR carousels rendered as literal escaped-tag text instead
  of an image slider.** Root cause: Pandoc's markdown reader treats any line
  indented ≥4 spaces as an indented code block, even inside an
  otherwise-open raw-HTML block — the old `carousel()` in `build_chapters.py`
  indented `<div class="carousel-item">`/`<img>` lines by 4-6 spaces. Also,
  inline elements (`<img>`, `<span>`) placed alone on their own line get
  wrapped in a stray `<p>` by pandoc. Fixed by flattening the whole carousel
  to zero indentation and keeping each inline element on the same line as its
  enclosing block tag. Verified with `pandoc -f markdown -t html` on
  Bsep/Cyp2c39/Cyp2c50 (multi-page carousels) — no `<pre><code>` escaping.
- Duplicate/stale figure trees (`figs/`, top-level `figures/`, a nested
  `figures/figures/` copy) merged into the single `figures/{protein}/` tree
  that `run_final_fig.sh` writes to.
- TIMBR sections were `include=FALSE`/commented out in every chapter despite
  real output existing for 26/27 proteins (all but Ces2) — now wired in.
