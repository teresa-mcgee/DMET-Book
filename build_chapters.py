#!/usr/bin/env python3
"""Deterministic generator for the DMET-results protein chapters and status.qmd.

Replaces the old GPT-driven generate_quartos.py as the source of truth for chapter
*structure*. It keeps each protein's existing "Protein Overview" paragraph (already
authored/reviewed) but rebuilds every other section by checking, on disk, which
figures actually exist for that protein -- so a chapter never shows a broken image
or claims an analysis is done when it isn't.

Run this after regenerating figures, or after re-running an analysis stage for a
protein, to refresh that protein's chapter and the Analysis Status table.

Usage: python3 build_chapters.py
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent          # DMET-results/
DMETPAPER = ROOT.parent                          # outer project root

PROTEINS = [
    l.strip().strip('"')
    for l in (ROOT / "biomarkernames_pyth.txt").read_text().splitlines()
    if l.strip()
]

SCAN_PAGE1_CANDIDATES = [
    "scans_0.95_page_001.jpg",
    "scans_0.9_page_001.jpg",
    "bc_scans_0.9_page_001.jpg",
    "scans_page_001.jpg",
    "allchScanV2_page_001.jpg",
]


def protein_dir(protein):
    return ROOT / "figures" / protein


def build_filename_index(pdir):
    """Walk pdir/jpgs ONCE and index every file by its basename. NFS-friendly:
    a naive rglob() per filename re-walks the whole subtree each call, which is
    far too slow across proteins like Cyp2c39/Cyp2c50 with thousands of TIMBR
    circos images."""
    jpgs = pdir / "jpgs"
    index = {}
    if not jpgs.is_dir():
        return index
    import os
    for dirpath, _dirnames, filenames in os.walk(jpgs):
        for fn in filenames:
            index.setdefault(fn, []).append(Path(dirpath) / fn)
    for fn in index:
        index[fn].sort()
    return index


def find_first(index, *filenames):
    """First existing file (by sorted path) matching any candidate filename, in priority order."""
    for name in filenames:
        matches = index.get(name)
        if matches:
            return matches[0]
    return None


def sibling_pages(page1_path):
    """All pages in the same series as a given `..._page_001.jpg` path."""
    if page1_path is None:
        return []
    m = re.match(r"^(.*)_page_\d+\.jpg$", page1_path.name)
    if not m:
        return [page1_path]
    stem = m.group(1)
    return sorted(
        page1_path.parent.glob(f"{stem}_page_*.jpg"),
        key=lambda p: int(re.search(r"_page_(\d+)\.jpg$", p.name).group(1)),
    )


def rel(path):
    return str(path.relative_to(ROOT))


def scan_variant_label(path):
    parts = {p.lower() for p in path.parts}
    if "boxcox" in parts:
        transform = "box-cox"
    elif "rint" in parts:
        transform = "RINT"
    elif "full" in parts:
        transform = "full dataset"
    else:
        transform = "transform unspecified"
    return f"{path.parent.name}, {transform}"


def count_scan_buckets(index):
    """Rough count of distinct scan variants (boxcox/rint/full) available, for the status table."""
    buckets = set()
    for name in SCAN_PAGE1_CANDIDATES:
        for m in index.get(name, []):
            parts = {p.lower() for p in m.parts}
            if "boxcox" in parts:
                buckets.add("boxcox")
            elif "rint" in parts:
                buckets.add("rint")
            elif "full" in parts:
                buckets.add("full")
            else:
                buckets.add("unlabeled")
    return len(buckets)


def rna_mediation_rds_exists(protein):
    candidates = [
        DMETPAPER / "output" / f"{protein}rna_mediators.RDS",
        DMETPAPER / "output" / "rna_mediators" / f"{protein}rna_mediators.RDS",
    ]
    return any(p.exists() for p in candidates)


def merge_status(protein):
    for base in [DMETPAPER / "merge" / protein, DMETPAPER / "out" / "Merge" / protein]:
        if base.is_dir():
            has_real_result = any(base.glob("mergeplot*")) or any(base.glob("*sigSNP_info*"))
            return "done" if has_real_result else "partial"
    return "not_run"


def extract_overview(protein):
    qmd = ROOT / f"{protein}.qmd"
    if not qmd.exists():
        return None
    text = qmd.read_text()
    m = re.search(r"# Protein Overview\s*\n+(.*?)\n#", text, re.DOTALL)
    return m.group(1).strip() if m else None


def img(path, width="80%"):
    return f'![]({rel(path)}){{fig-align="center" width="{width}"}}'


def carousel(ident, pages):
    if not pages:
        return ""
    if len(pages) == 1:
        return img(pages[0], width="100%")
    cid = re.sub(r"[^A-Za-z0-9_-]", "-", ident)
    items = []
    for i, p in enumerate(pages):
        active = " active" if i == 0 else ""
        items.append(f'<div class="carousel-item{active}"><img src="{rel(p)}" class="d-block w-100" alt="{cid} page {i + 1}"></div>')
    # Pandoc's markdown reader treats any line indented >=4 spaces as an indented
    # code block, even inside an otherwise-open raw-HTML block -- so this whole
    # block must stay at zero indentation, and each inline element (img, span)
    # must share a line with its enclosing block tag rather than sit on its own
    # line, or pandoc wraps it in a stray <p>. Previously indented, which rendered
    # as literal escaped tag text instead of an actual carousel.
    items_html = "\n".join(items)
    return f'''<div id="{cid}-carousel" class="carousel slide" data-bs-ride="carousel">
<div class="carousel-inner">
{items_html}
</div>
<button class="carousel-control-prev" type="button" data-bs-target="#{cid}-carousel" data-bs-slide="prev"><span class="carousel-control-prev-icon" aria-hidden="true"></span><span class="visually-hidden">Previous</span></button>
<button class="carousel-control-next" type="button" data-bs-target="#{cid}-carousel" data-bs-slide="next"><span class="carousel-control-next-icon" aria-hidden="true"></span><span class="visually-hidden">Next</span></button>
</div>'''


def missing(what, protein):
    return f"_Not yet analyzed here: no {what} found for {protein}. See [Analysis Status](status.qmd)._"


def analyze_protein(protein):
    """Everything the chapter builder AND the status table need, computed once."""
    pdir = protein_dir(protein)
    index = build_filename_index(pdir)
    raw = find_first(index, "rawHistograms_page_001.jpg")
    bc_hist = find_first(index, "BC_Histograms_page_001.jpg")
    rint_hist = find_first(index, "RINT_Histograms_page_001.jpg")
    boxcox_diag = find_first(index, "boxcox_page_001.jpg")
    h2 = pdir / "h2.jpg"
    combined = find_first(index, "genome_scan_combined_RNAProt.jpg", "genome_scan.jpg")
    scan_page1 = find_first(index, *SCAN_PAGE1_CANDIDATES)
    ci = find_first(index, "ci_scans_page_001.jpg")
    rna_med = find_first(index, "genome_wide_rna_med_page_001.jpg", "genome_wide_mediators_page_001.jpg")
    prot_med = find_first(index, "proteinMediators_page_001.jpg")
    timbr = find_first(index, "timbr_ci_page_001.jpg", "timbr_results_page_001.jpg")
    scan_buckets = count_scan_buckets(index)
    rna_rds = rna_mediation_rds_exists(protein)
    mstatus = merge_status(protein)

    n_dist = sum(1 for x in [raw, bc_hist, rint_hist, boxcox_diag] if x)
    distributions = "done" if n_dist == 4 else ("partial" if n_dist > 0 else "missing")
    heritability = "done" if h2.exists() else "missing"
    if scan_buckets >= 2 and ci:
        pqtl = "done"
    elif scan_buckets >= 1:
        pqtl = "partial"
    else:
        pqtl = "missing"
    if rna_med:
        rna = "done"
    elif rna_rds:
        rna = "partial"
    else:
        rna = "missing"
    timbr_status = "done" if timbr else "missing"

    return dict(
        pdir=pdir, raw=raw, bc_hist=bc_hist, rint_hist=rint_hist, boxcox_diag=boxcox_diag,
        h2=h2 if h2.exists() else None, combined=combined, scan_page1=scan_page1, ci=ci,
        rna_med=rna_med, prot_med=prot_med, timbr=timbr,
        distributions=distributions, heritability=heritability, pqtl=pqtl,
        rna_mediation=rna, merge=mstatus, timbr_status=timbr_status,
    )


def build_chapter(protein, a):
    overview = extract_overview(protein) or (
        f"_No protein overview text found on disk for {protein}; needs manual entry._"
    )

    L = []
    L += ['---', f'title: "{protein}"', 'format: html', 'editor: visual', '---', '']
    L += ['# Protein Overview', '', overview, '']
    L += ['# Data Overview', '',
          'The protein data can be handled either with RINT normalization or a Box-Cox (BC) '
          'transformation. BC transformation often produces clearer genome scans, but is not '
          'available yet for every analysis stage below.', '']

    L += ['## Protein Distributions', '']
    if a['boxcox_diag']:
        L += ['### Box-Cox transformation diagnostic', '', img(a['boxcox_diag']), '']
    for label, key in [("Raw data", "raw"), ("Box-Cox transformed", "bc_hist"), ("RINT transformed", "rint_hist")]:
        L += [f'### {label}', '']
        L += [img(a[key])] if a[key] else [missing(f"{label.lower()} distribution figure", protein)]
        L += ['']

    L += ['## Heritability', '',
          'Estimated narrow-sense heritability, i.e. whether the trait is strongly correlated '
          'within Collaborative Cross strains.', '']
    if a['h2']:
        L += [img(a['h2'])]
    else:
        L += [missing("heritability figure", protein)]
    L += ['']

    L += ['## pQTL Genome Scans', '']
    if a['combined']:
        L += ['### Combined with RNA (eQTL + pQTL)', '', img(a['combined']), '']
    if a['scan_page1']:
        pages = sibling_pages(a['scan_page1'])
        L += [f"### Genome-wide scan ({scan_variant_label(a['scan_page1'])})", '',
              carousel(protein, pages), '']
    else:
        L += [missing("pQTL genome scan figure", protein), '']

    L += ['## Peak Locus', '', '### Confidence interval of the significant peak (if one exists)', '']
    L += [img(a['ci'])] if a['ci'] else [missing("confidence-interval scan", protein)]
    L += ['']

    L += ['## RNA Mediation', '',
          'Genome-wide test of whether transcript abundance mediates the protein QTL (bmediatR).', '']
    if a['rna_med']:
        L += [img(a['rna_med'])]
    elif a['rna_mediation'] == 'partial':
        L += [f'_Analysis has been run for {protein} but no genome-wide mediation figure was '
              f'produced yet -- see `output/{protein}rna_mediators.RDS` (or `output/rna_mediators/`) '
              f'in the main project repository._']
    else:
        L += [missing("RNA-mediation figure", protein)]
    L += ['', '### Protein-protein mediation at the peak locus', '']
    L += [img(a['prot_med'])] if a['prot_med'] else [missing("protein-protein mediation figure", protein)]
    L += ['']

    L += ['## TIMBR (Allelic Series)', '',
          'TIMBR infers whether founder haplotypes at a QTL collapse into fewer functional '
          'alleles, using a Chinese Restaurant Process prior over haplotype groupings.', '']
    if a['timbr']:
        pages = sibling_pages(a['timbr'])
        L += [carousel(f"{protein}-timbr", pages) if len(pages) > 1 else img(a['timbr'])]
    else:
        L += [missing("TIMBR output", protein)]
    L += ['']

    L += ['## Merge / Diplotype Fine-Mapping', '']
    if a['merge'] == 'done':
        L += [f'Merge analysis has been run for {protein}; see `merge/{protein}/` in the main '
              f'project repository for diplotype and SNP-filtering results.']
    elif a['merge'] == 'partial':
        L += [f'A partial merge run exists for {protein} (diplotype filtering only, no '
              f'significant SNPs surfaced yet); see `out/Merge/{protein}/`.']
    else:
        L += [f'_Not yet analyzed here: merge/diplotype fine-mapping has not been run for '
              f'{protein}. See [Analysis Status](status.qmd)._']
    L += ['']

    return "\n".join(L)


STATUS_SYMBOL = {
    "done": "✅", "placeholder": "⚠️", "partial": "⚠️", "missing": "❌", "not_run": "❌",
}


def build_status_table(analyses):
    header = "| Protein | Distributions | Heritability* | pQTL scan | RNA mediation | Merge/diplotype | TIMBR |"
    sep = "|---|---|---|---|---|---|---|"
    rows = [header, sep]
    for protein in PROTEINS:
        a = analyses[protein]
        rows.append(
            f"| {protein} "
            f"| {STATUS_SYMBOL[a['distributions']]} "
            f"| {STATUS_SYMBOL[a['heritability']]} "
            f"| {STATUS_SYMBOL[a['pqtl']]} "
            f"| {STATUS_SYMBOL[a['rna_mediation']]} "
            f"| {STATUS_SYMBOL[a['merge']]} "
            f"| {STATUS_SYMBOL[a['timbr_status']]} |"
        )
    return "\n".join(rows)


def main():
    analyses = {p: analyze_protein(p) for p in PROTEINS}

    for protein, a in analyses.items():
        chapter = build_chapter(protein, a)
        (ROOT / f"{protein}.qmd").write_text(chapter + "\n")
        print(f"wrote {protein}.qmd")

    table = build_status_table(analyses)
    n_done = lambda key, val: sum(1 for a in analyses.values() if a[key] == val)
    print("\n--- summary ---")
    print(f"distributions done: {n_done('distributions', 'done')}/27")
    print(f"pQTL done: {n_done('pqtl', 'done')}/27, partial: {n_done('pqtl', 'partial')}/27")
    print(f"RNA mediation done: {n_done('rna_mediation', 'done')}/27, partial: {n_done('rna_mediation', 'partial')}/27")
    print(f"merge done: {n_done('merge', 'done')}/27, partial: {n_done('merge', 'partial')}/27")
    print(f"TIMBR done: {n_done('timbr_status', 'done')}/27")

    status_path = ROOT / "status.qmd"
    start, end = "<!-- STATUS_TABLE_START -->", "<!-- STATUS_TABLE_END -->"
    text = status_path.read_text()
    if start in text and end in text:
        i, j = text.index(start), text.index(end) + len(end)
        status_path.write_text(f"{text[:i]}{start}\n\n{table}\n\n{end}{text[j:]}")
        print("\nupdated status.qmd's table in place")
    else:
        print(f"\nWARNING: status.qmd is missing {start}/{end} markers -- table not written. "
              "Table is below to paste in manually:\n")
        print(table)


if __name__ == "__main__":
    main()
