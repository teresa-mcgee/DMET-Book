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

# Genome-wide scan carousels are page 1 = whole-genome overview, then one page per
# chromosome in this order (verified against Bsep/Cyp2c39/Ces2: page = chr index + 2).
CHR_ORDER = [str(i) for i in range(1, 20)] + ["X"]


def load_significant_peaks():
    """protein -> (chromosome, locus) of its most significant pQTL peak, from the
    permutation/GEV-adjusted genome-wide table (out/SignificantResults.txt) --
    NOT the per-protein output/{protein}/significantResults_continous.txt files,
    which use an uncorrected p<=0.05 cutoff and are too liberal (e.g. 10 "hits"
    for Oatp1a4) to trust for this."""
    path = DMETPAPER / "out" / "SignificantResults.txt"
    best = {}
    if not path.exists():
        return best
    lines = path.read_text().splitlines()[1:]
    for line in lines:
        parts = line.split()
        if len(parts) < 5:
            continue
        protein, chrom, locus, pval = parts[0], parts[1], parts[2], parts[4]
        try:
            pval = float(pval)
        except ValueError:
            continue
        if protein not in best or pval < best[protein][2]:
            best[protein] = (chrom, locus, pval)
    return {p: (chrom, locus) for p, (chrom, locus, _) in best.items()}


def load_all_peaks():
    """protein -> ALL its significant peaks (chrom, locus, pval), sorted by
    pval ascending (most significant first) -- unlike load_significant_peaks()
    above, which keeps only the single best one. Used for proteins with more
    than one significant locus on different chromosomes (currently only
    Ces2: chr18 and chr3), to order/label their per-chromosome TIMBR
    sections by significance rather than an arbitrary chromosome-number sort."""
    path = DMETPAPER / "out" / "SignificantResults.txt"
    all_peaks = {}
    if not path.exists():
        return all_peaks
    for line in path.read_text().splitlines()[1:]:
        parts = line.split()
        if len(parts) < 5:
            continue
        protein, chrom, locus, pval = parts[0], parts[1], parts[2], parts[4]
        try:
            pval = float(pval)
        except ValueError:
            continue
        all_peaks.setdefault(protein, []).append((chrom, locus, pval))
    for protein in all_peaks:
        all_peaks[protein].sort(key=lambda t: t[2])
    return all_peaks


SIGNIFICANT_PEAKS = load_significant_peaks()
SIGNIFICANT_CHR = {p: chrom for p, (chrom, _locus) in SIGNIFICANT_PEAKS.items()}
SIGNIFICANT_LOCUS = {p: locus for p, (_chrom, locus) in SIGNIFICANT_PEAKS.items()}
ALL_PEAKS = load_all_peaks()


def scan_active_index(protein, n_pages):
    """0-based index into the scan carousel's pages to open on, given this
    protein's significant chromosome. Falls back to 0 (first page) if unknown
    or out of range for however many pages this protein actually has."""
    chrom = SIGNIFICANT_CHR.get(protein)
    if chrom is None or chrom not in CHR_ORDER:
        return 0
    page_number = CHR_ORDER.index(chrom) + 2  # 1-based; page 1 is the overview
    index = page_number - 1
    return index if 0 <= index < n_pages else 0


# Transform subdirectories under out/{protein}/ that might hold a bootstrap-CI
# loci file, and the filename-separator variants seen historically (some
# scripts wrote "{locus}_ci.txt", others "{locus}ci.txt" or "{locus}_bc_ci.txt").
CI_FILE_DIRS = ["full", "boxcox", "rint", "bin", "sub"]


def find_ci_file(protein, locus):
    """Locate the bootstrap-CI loci file for a given protein/peak locus."""
    for d in CI_FILE_DIRS:
        base = DMETPAPER / "out" / protein / d
        for name in (f"{locus}_ci.txt", f"{locus}ci.txt", f"{locus}_bc_ci.txt"):
            candidate = base / name
            if candidate.exists() and candidate.stat().st_size > 0:
                return candidate
    return None


def timbr_peak_page_index(protein, locus, n_pages):
    """1-based page index of the peak locus within a `timbr_ci_page_*` carousel
    (one page per CI-sweep locus), used to power the TIMBR "Lead SNP"/5'/3'
    quick-jump buttons -- or None if it can't be determined with confidence.

    Deliberately conservative: only returns an index when the CI file's row
    count matches n_pages exactly and the peak locus is found in it, so a
    mismatched/mislabeled CI file (as happened historically for Cyp2c39 --
    see R/TIMBR/run_timbr_significant_gaps.R header) fails silently (no
    quick-jump buttons shown) rather than pointing at the wrong page."""
    ci_path = find_ci_file(protein, locus)
    if ci_path is None:
        return None
    rows = []
    for line in ci_path.read_text().splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            mb = float(parts[0])
        except ValueError:
            continue
        rows.append((mb, parts[2]))
    if len(rows) != n_pages:
        return None
    rows.sort(key=lambda r: r[0])  # 5' -> 3' along the chromosome
    for i, (_mb, loc) in enumerate(rows):
        if loc == locus:
            return i + 1  # 1-based page number
    return None


# Raw TIMBR sweep/peak images, as written directly by the CI-sweep loop in
# R/dev/final_figure_generation.R (~lines 358-450) and, for Ces2's
# single-locus fallback, R/TIMBR/run_timbr_significant_gaps.R
# (plot_singlelocus_peak()): under jpgs/{transform}/CHR{chrom}/TIMBR/ there is
# always effects/peak.jpg + circos/peak.jpg (the headline peak-locus result),
# and OPTIONALLY effects/index_N_{locus}.jpg + circos/index_N_{locus}.jpg for
# N = 1..len(sweep) (one pair per CI-sweep locus, absent for Ces2's
# single-locus fallback). This is a *different* on-disk convention from the
# `timbr_ci_page_NNN.jpg` composited-page series (Cyp2c39/Cyp2c50, made by a
# separate one-off compositing step) -- both exist for Cyp2c39/Cyp2c50, and
# the composite is deliberately preferred for those (see build_chapter()).
TIMBR_INDEX_RE = re.compile(r"^index_(\d+)_(.+)\.jpg$")


def find_timbr_raw(protein, pdir):
    """All `CHR{chrom}/TIMBR` raw-image dirs for this protein, across every
    transform subdirectory, as a {chrom: {...}} dict. Each value has
    'peak_effects'/'peak_circos' (Path or None) and 'sweep' (list of
    (index, locus, effects_path, circos_path) tuples, sorted by index
    ascending -- verified for Bsep/Ent1 that ascending index order equals
    ascending genomic (Mb) order, by cross-checking against the
    out/{protein}/{transform}/chr_{chrom}_ci_interval.txt loci list; empty
    list if this chrom has no sweep, just a single peak (Ces2)).

    If more than one transform has a TIMBR dir for the same chromosome
    (hasn't happened as of 2026-08), the one with more sweep loci wins, on
    the theory that a fuller sweep is more informative than a thinner one."""
    jpgs = pdir / "jpgs"
    if not jpgs.is_dir():
        return {}
    by_chrom = {}
    for timbr_dir in jpgs.glob("*/CHR*/TIMBR"):
        chrom = timbr_dir.parent.name[len("CHR"):]
        effects_dir, circos_dir = timbr_dir / "effects", timbr_dir / "circos"
        if not (effects_dir.is_dir() and circos_dir.is_dir()):
            continue
        sweep = []
        for f in effects_dir.glob("index_*_*.jpg"):
            m = TIMBR_INDEX_RE.match(f.name)
            if not m:
                continue
            circos_f = circos_dir / f.name
            if circos_f.exists():
                sweep.append((int(m.group(1)), m.group(2), f, circos_f))
        sweep.sort(key=lambda t: t[0])
        peak_effects = effects_dir / "peak.jpg"
        peak_circos = circos_dir / "peak.jpg"
        entry = dict(
            peak_effects=peak_effects if peak_effects.exists() else None,
            peak_circos=peak_circos if peak_circos.exists() else None,
            sweep=sweep,
        )
        if chrom not in by_chrom or len(sweep) > len(by_chrom[chrom]["sweep"]):
            by_chrom[chrom] = entry
    return by_chrom


def img_pair(path1, path2, width="48%"):
    """Two images side by side (haplotype-effects + circos plot), wrapped so
    they stack on narrow viewports. Zero-indentation raw HTML for the same
    pandoc reason documented in carousel()."""
    return (f'<div class="d-flex flex-wrap justify-content-center gap-2">'
            f'<img src="{rel(path1)}" style="max-width:{width};height:auto;">'
            f'<img src="{rel(path2)}" style="max-width:{width};height:auto;"></div>')


def timbr_sweep_carousel(ident, sweep, peak_locus=None):
    """Carousel over a CI-sweep's (index, locus, effects_path, circos_path)
    tuples, one page per locus showing both images side by side, with
    5'/Lead SNP/3' quick-jump buttons if peak_locus is found in the sweep."""
    pages = [(eff, cir) for _i, _loc, eff, cir in sweep]
    labels = [f"Figure {i + 1} ({loc})" for i, (_idx, loc, _e, _c) in enumerate(sweep)]
    active_index = 0
    quick_jumps = None
    if peak_locus is not None:
        for i, (_idx, loc, _e, _c) in enumerate(sweep):
            if loc == peak_locus:
                active_index = i
                quick_jumps = [("5' side", 0), ("Lead SNP", i), ("3' side", len(sweep) - 1)]
                break
    return carousel(ident, pages, active_index=active_index, quick_jumps=quick_jumps, labels=labels)


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


def carousel(ident, pages, active_index=0, quick_jumps=None, labels=None):
    """quick_jumps: optional list of (label, 0-based page index) tuples,
    rendered as buttons that jump straight to that page (e.g. TIMBR's
    5' side / Lead SNP / 3' side) -- in addition to the always-present
    "jump to figure" dropdown, which covers the general "jump to the Nth
    figure" case for every multi-page carousel.

    Each element of `pages` is normally a single image Path, but may instead
    be a (path1, path2) tuple -- e.g. a TIMBR sweep's haplotype-effects +
    circos-plot pair for one locus -- rendered side by side within that page.

    labels: optional list of per-page strings for the "jump to figure"
    dropdown (defaults to "Figure N")."""
    if not pages:
        return ""
    if len(pages) == 1:
        p = pages[0]
        return img_pair(*p) if isinstance(p, tuple) else img(p, width="100%")
    if not (0 <= active_index < len(pages)):
        active_index = 0
    cid = re.sub(r"[^A-Za-z0-9_-]", "-", ident)
    items = []
    for i, p in enumerate(pages):
        active = " active" if i == active_index else ""
        if isinstance(p, tuple):
            inner = "".join(
                f'<img src="{rel(x)}" style="max-width:48%;height:auto;">' for x in p
            )
            items.append(f'<div class="carousel-item{active}"><div class="d-flex flex-wrap justify-content-center gap-2">{inner}</div></div>')
        else:
            items.append(f'<div class="carousel-item{active}"><img src="{rel(p)}" class="d-block w-100" alt="{cid} page {i + 1}"></div>')
    # Pandoc's markdown reader treats any line indented >=4 spaces as an indented
    # code block, even inside an otherwise-open raw-HTML block -- so this whole
    # block must stay at zero indentation, and each inline element (img, span)
    # must share a line with its enclosing block tag rather than sit on its own
    # line, or pandoc wraps it in a stray <p>. Previously indented, which rendered
    # as literal escaped tag text instead of an actual carousel. The jump
    # controls below follow the same rule.
    items_html = "\n".join(items)

    quick_jump_html = ""
    if quick_jumps:
        buttons = [
            f'<button type="button" class="btn btn-sm btn-outline-secondary" data-bs-target="#{cid}-carousel" data-bs-slide-to="{idx}">{label}</button>'
            for label, idx in quick_jumps
        ]
        quick_jump_html = "\n".join(buttons) + "\n"

    options = []
    for i in range(len(pages)):
        selected = " selected" if i == active_index else ""
        label = labels[i] if labels else f"Figure {i + 1}"
        options.append(f'<option value="{i}"{selected}>{label}</option>')
    options_html = "\n".join(options)

    jump_controls = f'''<div class="carousel-jump-controls d-flex flex-wrap align-items-center gap-2 mt-2">
{quick_jump_html}<label for="{cid}-jump" class="mb-0">Jump to figure:</label>
<select id="{cid}-jump" class="form-select form-select-sm d-inline-block w-auto" onchange='bootstrap.Carousel.getOrCreateInstance(document.getElementById("{cid}-carousel")).to(parseInt(this.value));'>
{options_html}
</select>
</div>
<script>document.getElementById("{cid}-carousel").addEventListener("slide.bs.carousel",function(e){{document.getElementById("{cid}-jump").value=e.to;}});</script>'''

    return f'''<div id="{cid}-carousel" class="carousel slide" data-bs-ride="carousel">
<div class="carousel-inner">
{items_html}
</div>
<button class="carousel-control-prev" type="button" data-bs-target="#{cid}-carousel" data-bs-slide="prev"><span class="carousel-control-prev-icon" aria-hidden="true"></span><span class="visually-hidden">Previous</span></button>
<button class="carousel-control-next" type="button" data-bs-target="#{cid}-carousel" data-bs-slide="next"><span class="carousel-control-next-icon" aria-hidden="true"></span><span class="visually-hidden">Next</span></button>
</div>
{jump_controls}'''


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
    timbr_raw = find_timbr_raw(protein, pdir)
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
    has_raw_peak = any(e["peak_effects"] or e["peak_circos"] for e in timbr_raw.values())
    timbr_status = "done" if (timbr or has_raw_peak) else "missing"

    return dict(
        pdir=pdir, raw=raw, bc_hist=bc_hist, rint_hist=rint_hist, boxcox_diag=boxcox_diag,
        h2=h2 if h2.exists() else None, combined=combined, scan_page1=scan_page1, ci=ci,
        rna_med=rna_med, prot_med=prot_med, timbr=timbr, timbr_raw=timbr_raw,
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
        active_index = scan_active_index(protein, len(pages))
        L += [f"### Genome-wide scan ({scan_variant_label(a['scan_page1'])})", '']
        if active_index:
            L += [f'_Opens on Chr {SIGNIFICANT_CHR[protein]}, where the significant peak is._', '']
        L += [carousel(protein, pages, active_index=active_index), '']
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
    L += [img(a['prot_med'], width="100%")] if a['prot_med'] else [missing("protein-protein mediation figure", protein)]
    L += ['']

    L += ['## TIMBR (Allelic Series)', '',
          'TIMBR infers whether founder haplotypes at a QTL collapse into fewer functional '
          'alleles, using a Chinese Restaurant Process prior over haplotype groupings.', '']
    timbr_raw = a['timbr_raw']
    # The raw per-locus images (effects/peak.jpg + circos/peak.jpg, optionally
    # with an index_N_{locus} sweep) are preferred over a `timbr_results_page_*`
    # composite (a stale/mismatched 2-page RINT summary, currently Bsep/Ent1)
    # -- but NOT over a `timbr_ci_page_*` composite (Cyp2c39/Cyp2c50), which is
    # already a faithful, working rendering of the same underlying sweep and
    # is left untouched.
    use_raw = bool(timbr_raw) and not (a['timbr'] and a['timbr'].name.startswith('timbr_ci_page'))
    if use_raw:
        # Order multi-peak proteins (currently only Ces2: chr18 + chr3) by
        # significance (most significant first); single-peak proteins have
        # only one chrom here so this is a no-op for them.
        ranked = [c for c, _l, _p in ALL_PEAKS.get(protein, []) if c in timbr_raw]
        chroms = ranked + [c for c in sorted(timbr_raw) if c not in ranked]
        multi = len(chroms) > 1
        for chrom in chroms:
            entry = timbr_raw[chrom]
            if multi:
                L += [f'### Chr {chrom} peak locus', '']
            if entry['peak_effects'] and entry['peak_circos']:
                L += [img_pair(entry['peak_effects'], entry['peak_circos']), '']
            elif entry['peak_effects'] or entry['peak_circos']:
                L += [img(entry['peak_effects'] or entry['peak_circos'], width="60%"), '']
            else:
                L += [missing(f"TIMBR peak-locus figure for Chr {chrom}", protein), '']
            if entry['sweep']:
                # Quick-jump buttons only make sense when this chromosome's
                # peak matches the protein's single tracked significant locus
                # (SIGNIFICANT_LOCUS) -- for a multi-peak protein like Ces2,
                # neither chrom's peak is individually tracked there, so the
                # sweep (if one ever exists) would just open on page 1.
                peak_locus = SIGNIFICANT_LOCUS.get(protein) if not multi else None
                L += ["Haplotype effects and circos plot swept across the bootstrap confidence "
                      "interval (5' to 3' along the chromosome):", '']
                L += [timbr_sweep_carousel(f"{protein}-timbr-chr{chrom}", entry['sweep'], peak_locus), '']
            else:
                L += [f'_Single-locus TIMBR result only -- the bootstrap-CI sweep failed/was empty for '
                      f'this locus, so only the peak marker itself was tested. See '
                      f'[Analysis Status](status.qmd)._', '']
    elif a['timbr']:
        pages = sibling_pages(a['timbr'])
        quick_jumps = None
        # Only the true per-CI-locus sweep series (`timbr_ci_page_*`, one page
        # per bootstrap-CI locus, currently Cyp2c39/Cyp2c50) has a meaningful
        # 5'/3'/lead-SNP structure -- `timbr_results_page_*` is a generic
        # 2-page haplotype+circos summary of the single peak locus, so no
        # quick-jump buttons are offered there.
        if len(pages) > 1 and a['timbr'].name.startswith("timbr_ci_page") and protein in SIGNIFICANT_LOCUS:
            peak_idx1 = timbr_peak_page_index(protein, SIGNIFICANT_LOCUS[protein], len(pages))
            if peak_idx1:
                quick_jumps = [
                    ("5' side", 0),
                    ("Lead SNP", peak_idx1 - 1),
                    ("3' side", len(pages) - 1),
                ]
        L += [carousel(f"{protein}-timbr", pages, quick_jumps=quick_jumps) if len(pages) > 1 else img(a['timbr'])]
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
