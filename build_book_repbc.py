#!/usr/bin/env python3
"""Rebuild the whole DMET-results book from the replicate-level Box-Cox analysis (tag repbc_2026-10-01).

Generates: index.qmd, intro.qmd (Methods), overview.qmd, status.qmd, credible_intervals.qmd and the 26
protein chapters. Every figure is shown INLINE (not just linked) and followed by a "Download" line that
offers the PNG / PDF (or CSV) as a file download. The hand-written "Protein Overview" paragraph of each
chapter is preserved. Anything that has not been produced yet (e.g. Merge / TIMBR plots that are still
running) is reported as such rather than silently omitted -- re-run this script after those finish.

Usage (from DMET-results/):  python3 build_book_repbc.py
Inputs: ../out/<TAG>/, ../out/<TAG>_zeroinflated/, ../merge/<TAG>/, figures/<TAG>/ (made by R/rerun2026/*.R).
After running, render with:  quarto render
"""
import csv
import datetime
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
TAG = "repbc_2026-10-01"
ZTAG = TAG + "_zeroinflated"
RNA_TAG = "rerun_2026-09-27_rna_v2"
OUT, ZOUT = PROJECT / "out" / TAG, PROJECT / "out" / ZTAG
MERGE, ZMERGE = PROJECT / "merge" / TAG, PROJECT / "merge" / ZTAG
FIGS = ROOT / "figures" / TAG
PP = FIGS / "per_protein"
PROTEINS = [l.strip().strip('"') for l in (ROOT / "biomarkernames_pyth.txt").read_text().splitlines() if l.strip()]
PROTEINS = [p for p in PROTEINS if p != "Ces2"]
ZERO_INFLATED = ["Cyp2c39", "Cyp2c50"]          # proteins with a zero-inflated subset scan (Ces2 excluded)
FIGS.mkdir(parents=True, exist_ok=True)
PP.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------- helpers
def load_csv(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def load_sig(path):
    rows = []
    for line in Path(path).read_text().splitlines()[1:]:
        t = line.split()
        if len(t) >= 10:
            rows.append(dict(Protein=t[0], CH=t[1], Locus=t[2], Mb=float(t[3]), P=float(t[4]), CI_start=float(t[5]),
                             CI_end=float(t[6]), Subset=t[7], Adj_P=float(t[9])))
    return [r for r in rows if r["Protein"] != "Ces2"]


def num(x):
    try:
        v = float(x)
        return None if v != v else v
    except (TypeError, ValueError):
        return None


def f2(x, d=2):
    v = num(x)
    return "—" if v is None else f"{v:.{d}f}"


def fp(x):
    v = num(x)
    return "—" if v is None else f"{v:.2e}"


def rel(p):
    return str(Path(p).relative_to(ROOT))


def stage(src, protein, name=None):
    """Copy a source file into the book tree (figures/<TAG>/per_protein/<Protein>/) and return the new path."""
    src = Path(src)
    fname = name or src.name
    if not fname.lower().startswith(protein.lower()):
        fname = f"{protein}_{fname}"
    dst = PP / protein / fname
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists() or dst.stat().st_mtime < src.stat().st_mtime:
        shutil.copy2(src, dst)
    return dst


def pdf_to_png(pdf, png, page=1, dpi=130):
    pdf, png = Path(pdf), Path(png)
    if png.exists() and png.stat().st_mtime >= pdf.stat().st_mtime:
        return png
    png.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["pdftoppm", "-png", "-r", str(dpi), "-f", str(page), "-l", str(page), "-singlefile", str(pdf),
                    str(png.with_suffix(""))], check=True)
    return png


def figure(png, caption, downloads=(), width="95%"):
    """Inline figure + a download line (each entry of `downloads` is a path to an existing file)."""
    png = Path(png)
    dl = []
    seen = set()
    for d in downloads:
        d = Path(d)
        if d.exists() and d not in seen:
            seen.add(d)
            dl.append(f'[{d.suffix.lstrip(".").upper()}]({rel(d)}){{download="{d.name}"}}')
    block = [f'![{caption}]({rel(png)}){{fig-align="center" width="{width}"}}', ""]
    if dl:
        block += ["::: {.fig-downloads}", "Download: " + " · ".join(dl), ":::", ""]
    return "\n".join(block)


def tfig(stem, caption, width="95%", extra=()):
    """Figure from the book's figure folder: <stem>.png inline, <stem>.png/.pdf (+extra) downloadable."""
    png = FIGS / f"{stem}.png"
    if not png.exists():
        return f"_Figure not available yet: `{stem}`._\n"
    return figure(png, caption, [png, FIGS / f"{stem}.pdf", *[FIGS / e for e in extra]], width)


def dl_line(label, path):
    path = Path(path)
    return f'[{label}]({rel(path)}){{download="{path.name}"}}'


def short_symbol(sym):
    parts = [x.strip() for x in str(sym).split("///")]
    return parts[0] if len(parts) == 1 else f"{parts[0]} +{len(parts) - 1}"


def md(x):
    return str(x).replace("|", "\\|") if x not in (None, "") else "—"


# --------------------------------------------------------------------------- data
sig = load_sig(OUT / "SignificantResults.txt")
zsig = load_sig(ZOUT / "SignificantResults.txt") if (ZOUT / "SignificantResults.txt").exists() else []
for r in zsig:
    r["Subset"] = "sub"
norm = {r["Protein"].upper(): r for r in load_csv(FIGS / "protein_normalization_summary.csv")}
her = {r["Protein"].upper(): r for r in load_csv(FIGS / "figure2_heritability_protein_and_transcript.csv")}
pairs = {r["Protein"].upper(): r for r in load_csv(OUT / "protein_transcript_pairs.csv")}
cis = {}
for src in (OUT, ZOUT):
    f = src / "pqtl_cis_trans.csv"
    if f.exists():
        for r in load_csv(f):
            cis[(r["Protein"].upper(), r["Locus"])] = r
lead = {r["Protein"].upper(): r for r in load_csv(OUT / "prot_summary_table_with_heritability.csv")}


def sig_for(protein, rows=None):
    return [r for r in (rows if rows is not None else sig) if r["Protein"].upper() == protein.upper()]


def merge_state(protein, locus, subset="full"):
    base = (ZMERGE if subset == "sub" else MERGE) / protein / locus
    if (base / "sigSNP_info_v3.txt").exists() or (base / "diplo_res_filter.txt").exists():
        return "complete"
    return "running or not started" if base.exists() else "not started"


def fig_dir(subset):
    return (ZOUT if subset == "sub" else OUT) / "figs"


# --------------------------------------------------------------------------- allelic-series text, UniProt
FOUNDERS = ["A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ"]
SUMMARY_FILE = OUT / "pqtl_summary_table.csv"          # written by R/rerun2026/summarize_pqtl_analysis.R (after Merge)
SUMMARY = {}
if SUMMARY_FILE.exists():
    for _r in load_csv(SUMMARY_FILE):
        SUMMARY[(_r["Protein"].upper(), _r["Locus"])] = _r
UNIPROT_CACHE = FIGS / "uniprot_mediator_function_cache.json"
MEDIATION_LABEL_BF = 0.5          # same rule as the plot: label transcripts with best-model log10 BF >= 0.5


def allele_groups(code):
    """'0,0,0,1,0,0,0,1' -> {'0': [founders...], '1': [...]} (None if the code is not 8 values)."""
    if not code or code in ("NA", "—"):
        return None
    t = [x.strip() for x in str(code).replace(";", ",").split(",") if x.strip() != ""]
    if len(t) != 8:
        return None
    g = {}
    for founder, a in zip(FOUNDERS, t):
        g.setdefault(a, []).append(founder)
    return g


def allele_text(code, what):
    g = allele_groups(code)
    if g is None:
        return None
    parts = [f"**allele {a}**: {', '.join(v)}" for a, v in sorted(g.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    return f"{what} " + "; ".join(parts) + f". That is {len(g)} functional allele{'s' if len(g) > 1 else ''}: founders in the same group are inferred to carry the same functional allele at this locus."


def merge_text(protein, locus):
    r = SUMMARY.get((protein.upper(), locus))
    if r is not None and r.get("merge_source") != "tag":      # SDP not read from THIS tag's Merge output -> never show it
        r = None
    if r is None:
        return "_The top allelic series (SDP) from Merge is added to the book once the follow-up summary has been generated._\n"
    t = allele_text(r.get("top_sdp"), "In plain terms, the strongest Merge strain distribution pattern (SDP) splits the eight founders into")
    return (t + "\n") if t else "_Merge did not return a top SDP for this locus (not run or no significant SDP)._\n"


def timbr_text(protein, locus):
    r = SUMMARY.get((protein.upper(), locus))
    if r is None:
        return None
    t = allele_text(r.get("top_timbr_allelic_series"), "The highest-posterior TIMBR allelic series at the peak splits the founders into")
    return (t + "\n") if t else None


def _clean_function(txt, limit=380):
    txt = re.sub(r"\s*\((?:PubMed|By similarity)[^)]*\)", "", txt or "")
    txt = re.sub(r"\s+", " ", txt).strip()
    if len(txt) <= limit:
        return txt
    cut = txt[:limit]
    k = max(cut.rfind(". "), cut.rfind("; "))
    return (cut[:k + 1] if k > 150 else cut.rsplit(" ", 1)[0] + " …").strip()


def uniprot_lookup(symbols):
    """Function annotation for mouse gene symbols from the UniProt REST API (cached on disk; reviewed entries first)."""
    cache = json.loads(UNIPROT_CACHE.read_text()) if UNIPROT_CACHE.exists() else {}
    todo = sorted({x for x in symbols if x and x not in cache})
    def fetch(batch, reviewed):
        q = "(" + " OR ".join(f"gene_exact:{x}" for x in batch) + f") AND organism_id:10090 AND reviewed:{'true' if reviewed else 'false'}"
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
            {"query": q, "fields": "accession,gene_names,protein_name,cc_function", "format": "json", "size": 500})
        with urllib.request.urlopen(url, timeout=60) as resp:
            return json.load(resp).get("results", [])
    def parse(e):
        names = [g.get("geneName", {}).get("value") for g in e.get("genes", [])] + [s.get("value") for g in e.get("genes", []) for s in g.get("synonyms", [])]
        pd = e.get("proteinDescription", {})
        pname = (pd.get("recommendedName") or (pd.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value", "")
        fn = [c["texts"][0]["value"] for c in e.get("comments", []) if c.get("commentType") == "FUNCTION" and c.get("texts")]
        return [x for x in names if x], dict(accession=e.get("primaryAccession"), protein=pname, function=_clean_function(fn[0]) if fn else "", reviewed=None)
    try:
        for reviewed in (True, False):
            pending = [x for x in todo if x not in cache or cache[x].get("status") != "found"]
            for i in range(0, len(pending), 40):
                batch = pending[i:i + 40]
                for e in fetch(batch, reviewed):
                    names, rec = parse(e)
                    rec["reviewed"] = reviewed
                    for nm in names:
                        if nm in batch and (nm not in cache or cache[nm].get("status") != "found" or (reviewed and not cache[nm].get("reviewed"))):
                            cache[nm] = dict(rec, status="found")
                time.sleep(0.3)
        for x in todo:
            cache.setdefault(x, dict(status="none"))
        UNIPROT_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
    except (urllib.error.URLError, TimeoutError, OSError) as err:     # offline build: use whatever is cached
        print(f"[uniprot] lookup skipped ({err}); using cached entries only")
    return cache


def mediation_points(protein, chr_, fdir):
    f = fdir / f"mediation_bf_points_chr{chr_}.csv"
    if not f.exists():
        return None, f
    rows = load_csv(f)
    return rows, f


def mediation_text_and_table(protein, chr_, fdir):
    rows, f = mediation_points(protein, chr_, fdir)
    if rows is None:
        return "", ""
    best = {}
    for r in rows:                                                    # one row per gene symbol: strongest RNA column
        b = float(r["best_log10BF"])
        if r["symbol"] not in best or b > float(best[r["symbol"]]["best_log10BF"]):
            best[r["symbol"]] = r
    lab = sorted([r for r in best.values() if r.get("labelled") in ("TRUE", "True", "true")],
                 key=lambda r: (r["best_model"] != "mediation", -float(r["best_log10BF"])))
    n_other = sum(1 for r in best.values() if r["best_model"] == "other" and float(r["best_log10BF"]) >= MEDIATION_LABEL_BF)
    names = {"mediation": "mediation", "colocal": "co-local", "other": "other (non-mediator)"}
    by = {m: [r for r in lab if r["best_model"] == m] for m in names}
    sent = (f"{len(best)} candidate transcripts lie in the credible interval. **{len(lab)}** have a best-model log10 Bayes factor ≥ {MEDIATION_LABEL_BF} with mediation or co-local as the best-supported model and are labelled in the figure"
            + (f" ({n_other} more reach {MEDIATION_LABEL_BF} only for the non-mediator model and are not labelled)" if n_other else "") + ": ")
    bits = []
    for m in ("mediation", "colocal"):
        if by[m]:
            bits.append(f"**{names[m]}** ({len(by[m])}): " + ", ".join(f"{r['symbol']} ({float(r['best_log10BF']):.2f})" for r in by[m][:12]) + (" …" if len(by[m]) > 12 else ""))
    sent += "; ".join(bits) if bits else "none."
    sent += ". *Mediation* means the transcript is best explained as lying on the path from the QTL to the protein; *co-local* means it shares the QTL without mediating; log10 BF is the evidence for the best-supported model against the null.\n"
    table = ""
    if lab:
        info = uniprot_lookup([r["symbol"] for r in lab])
        t = [f"::: {{.callout-note collapse=\"true\" title=\"Function of the {len(lab)} labelled transcripts (UniProt)\"}}",
             f"Function text is the UniProt (mouse) annotation retrieved on {datetime.date.today().isoformat()}; reviewed Swiss-Prot entries where available. "
             "A dash means UniProt holds no functional annotation for that symbol.\n",
             "| Transcript | Best model | log10 BF | UniProt protein | Function (UniProt) |", "|---|---|---:|---|---|"]
        for r in lab:
            u = info.get(r["symbol"], {})
            acc = u.get("accession")
            prot = f"[{md(u.get('protein'))}](https://www.uniprot.org/uniprotkb/{acc})" if acc and u.get("protein") else "—"
            t.append(f"| *{r['symbol']}* | {names[r['best_model']]} | {float(r['best_log10BF']):.2f} | {prot} | {md(u.get('function') or '')} |")
        t.append(":::\n")
        table = "\n".join(t)
    return sent, table


# --------------------------------------------------------------------------- per-locus figure blocks
def locus_blocks(protein, peak, level="##"):
    """Everything produced for one significant peak, shown inline."""
    subset = peak.get("Subset", "full")
    chr_, locus = peak["CH"], peak["Locus"]
    fdir = fig_dir(subset) / protein / ("sub" if subset == "sub" else "full")
    peaks = sig_for(protein, zsig if subset == "sub" else sig)
    j = peaks.index(peak) + 1
    ct = cis.get((protein.upper(), locus), {})
    out = []
    kind = ct.get("pqtl_type") or "—"
    gene = ct.get("pqtl_gene") or "—"
    gloc = f"chr{ct.get('gene_chr')}:{ct.get('gene_start_Mb')}–{ct.get('gene_end_Mb')} Mb" if ct.get("gene_chr") not in (None, "", "NA") else "—"
    out.append(f"{level} Chromosome {chr_} peak ({locus})" + (" — zero-inflated subset scan" if subset == "sub" else "") + "\n")
    out.append("| Peak position | −log10(P) | Adjusted P | 90% credible interval | Encoding gene | Gene location | pQTL type |\n|---|---:|---:|---|---|---|---|\n"
               f"| chr{chr_}:{peak['Mb']:.2f} Mb | {-math.log10(peak['P']):.2f} | {fp(peak['Adj_P'])} | "
               f"{peak['CI_start']:.2f}–{peak['CI_end']:.2f} Mb | {md(gene)} | {gloc} | **{kind}** |\n")
    if kind in ("cis", "trans"):
        out.append("A pQTL is **cis** when the gene encoding the protein lies within the credible interval and **trans** otherwise.\n")
    # bootstrap CI scan (page j of ci_scans.pdf) and gene-track CI plot
    cs = fdir / "ci_scans.pdf"
    if cs.exists():
        png = pdf_to_png(cs, PP / protein / f"{protein}_ci_scans_{subset}_p{j}.png", page=j)
        out.append(figure(png, f"Bootstrap credible interval, chromosome {chr_}.", [stage(cs, protein, f"ci_scans_{subset}.pdf")]))
    cg = fdir / f"ci_genes_chr{chr_}.png"
    if cg.exists():
        out.append(figure(stage(cg, protein, f"ci_genes_{subset}_chr{chr_}.png"),
                          f"Genes within the credible interval, chromosome {chr_} (80% interval boxed).", [stage(cg, protein, f"ci_genes_{subset}_chr{chr_}.png")]))
    # Merge
    mp = fdir / f"mergeplot_{locus}.png"
    out.append(f"{level}# Merge / diplotype fine-mapping\n")
    if mp.exists():
        out.append(figure(stage(mp, protein, f"mergeplot_{subset}_{locus}.png"), "Merge: SDP significance profile across the credible interval.",
                          [stage(mp, protein, f"mergeplot_{subset}_{locus}.png")]))
    else:
        out.append(f"_Merge status: {merge_state(protein, locus, subset)}; the plot is not available yet._\n")
    if subset == "full":
        out.append(merge_text(protein, locus))
    # TIMBR
    out.append(f"{level}# TIMBR allelic series\n")
    tdir = fig_dir(subset) / protein / "TIMBR_ci"
    timbr_imgs = sorted([p for p in tdir.glob("*") if p.suffix.lower() in {".png", ".jpg"} and f"chr{chr_}" in p.name and not p.name.startswith(("allele_number", "haplotype_"))]) if tdir.exists() else []
    if timbr_imgs:
        for p in timbr_imgs:
            kind = ("Bayes-factor" if "_bf_" in p.name else "posterior-probability")
            out.append(figure(stage(p, protein), f"TIMBR allelic-series sweep across the credible interval ({kind} dot plot), chromosome {chr_}; Merge SDP profile overlaid where available.", [stage(p, protein)]))
    else:
        done = (OUT / protein / "TIMBR_ci" / f"ch{chr_}_res_list.RData").exists() or (OUT / protein / "TIMBR_ci" / f"ch{chr_}_singlelocus_res.RData").exists()
        out.append("_TIMBR sweep complete; plots will appear once Merge profiles are available._\n" if (done and subset == "full")
                   else "_TIMBR was not run for the zero-inflated subset peak._\n" if subset == "sub" else "_TIMBR not run yet._\n")
    tt = timbr_text(protein, locus) if subset == "full" else None
    if tt:
        out.append(tt)
    if subset == "full":
        out.append(allele_number_blocks(protein, chr_, tdir))
    # mediation
    out.append(f"{level}# RNA mediation\n")
    if subset == "sub":
        out.append("_RNA mediation was not run for the zero-inflated subset peak._\n")
    else:
        bf_no = fdir / f"mediation_overlay_bf_chr{chr_}_nogenes.png"      # the book shows the plot WITHOUT the gene track
        bf_with = fdir / f"mediation_overlay_bf_chr{chr_}_withgenes.png"  # offered as a download
        pb = fdir / f"targeted_mediation_chr{chr_}_posterior_bars.pdf"
        pts_rows, pts_csv = mediation_points(protein, chr_, fdir)
        downloads = [stage(x, protein) for x in (bf_no, bf_with, pts_csv) if x.exists()]
        if bf_no.exists():
            out.append(figure(stage(bf_no, protein), "RNA mediation: pQTL scan (left axis) and each candidate transcript's best-supported model, log10 Bayes factor (right axis); transcripts with a mediation or co-local best model and log10 BF of at least 0.5 are labelled.", downloads))
            sent, table = mediation_text_and_table(protein, chr_, fdir)
            if sent:
                out.append(sent)
            if table:
                out.append(table)
        if pb.exists():
            png = pdf_to_png(pb, PP / protein / f"{protein}_mediation_posterior_bars_chr{chr_}.png", page=1)
            out.append(figure(png, "Posterior probabilities of the mediation models for the strongest candidate mediator.", [stage(pb, protein)]))
        if not bf_no.exists() and not pb.exists():
            out.append("_No mediation output yet for this peak._\n")
    return "\n".join(out)


ALLELE_SUMMARY = {}
_af = OUT / "timbr_allele_number_summary.csv"
if _af.exists():
    for _r in load_csv(_af):
        ALLELE_SUMMARY[(_r["Protein"].upper(), _r["CH"])] = _r


def allele_number_blocks(protein, chr_, tdir):
    """Number-of-functional-alleles (prior vs posterior) and lead / best-model haplotype plots for one peak."""
    r = ALLELE_SUMMARY.get((protein.upper(), str(chr_)))
    if not r:
        return ""
    def pk(pref, k):
        return float(r[f"{pref}_P_K{k}"])
    def ge4(pref):
        return sum(pk(pref, k) for k in range(4, 9))
    out = []
    out.append("**Number of functional alleles.** "
               f"At the lead locus ({r['lead_locus']}, {float(r['lead_Mb']):.2f} Mb) the posterior probability of a bi-allelic series (two functional alleles) is {pk('lead',2):.2f}, "
               f"tri-allelic {pk('lead',3):.2f}, and four or more alleles {ge4('lead'):.2f}; a single allele (no QTL effect) has probability {pk('lead',1):.2f}. "
               f"For comparison, the Chinese-restaurant-process prior gives {float(r['prior_crp_K1']):.2f} / {float(r['prior_crp_K2']):.2f} / {float(r['prior_crp_K3']):.2f} for one / two / three alleles, "
               f"so the data move weight away from a single allele. ")
    if r["lead_equals_best"] not in ("TRUE", "True", "1"):
        out[-1] += (f"The locus in the credible interval whose single best allelic series has the highest posterior probability is {r['best_locus']} ({float(r['best_Mb']):.2f} Mb; "
                    f"series {r['best_top_partition']}, P = {float(r['best_top_partition_P']):.2f}); there bi-allelic = {pk('best',2):.2f}, tri-allelic = {pk('best',3):.2f}, four or more = {ge4('best'):.2f}.")
    else:
        out[-1] += "The lead locus is also the locus with the highest-probability allelic series."
    out[-1] += "\n"
    items = [
        (f"allele_number_prior_posterior_chr{chr_}", "Prior (uniform over partitions and Chinese restaurant process) versus posterior probability of the number of functional alleles at the lead and best-model loci."),
        (f"allele_number_across_ci_chr{chr_}", "Posterior probability of 1, 2, 3 or 4+ functional alleles at each locus across the credible interval (dashed line: lead locus; dotted line: best-model locus)."),
        (f"haplotype_lead_{r['lead_locus']}_chr{chr_}", f"TIMBR founder-haplotype effects at the lead locus {r['lead_locus']}."),
        (f"haplotype_best_{r['best_locus']}_chr{chr_}", f"TIMBR founder-haplotype effects at the locus with the highest-probability allelic series, {r['best_locus']}."),
    ]
    for stem, cap in items:
        png = tdir / f"{stem}.png"
        if png.exists():
            out.append(figure(stage(png, protein), cap, [stage(png, protein), tdir / f"{stem}.pdf" if not (PP / protein / f"{protein}_{stem}.pdf").exists() else None] and [stage(png, protein), stage(tdir / f"{stem}.pdf", protein)] if (tdir / f"{stem}.pdf").exists() else [stage(png, protein)]))
    return "\n".join(out)


# --------------------------------------------------------------------------- protein chapters
def extract_overview(protein):
    f = ROOT / f"{protein}.qmd"
    if f.exists():
        m = re.search(r"# Protein Overview\n(.*?)\n# ", f.read_text(), re.S)
        if m:
            return m.group(1).strip()
    return "_No overview paragraph recorded for this protein yet._"


def build_chapter(protein):
    key = protein.upper()
    nr, hr, pr, ld = norm[key], her[key], pairs[key], lead.get(key, {})
    peaks = sig_for(protein)
    L = [f'---\ntitle: "{protein}"\nformat: html\neditor: visual\n---\n', "# Protein Overview\n", extract_overview(protein) + "\n"]
    # --- data
    L.append("# Data and normalization\n")
    shift = num(nr["shift"]) or 0
    L.append(f"Replicate measurements were Box-Cox transformed with a protein-specific λ = **{f2(nr['lambda'])}** "
             f"(95% likelihood interval {f2(nr['lambda_CI_low'])} to {f2(nr['lambda_CI_high'])}), estimated on the individual replicates "
             f"(`boxcox(y ~ strain)`, no +1 shift)" + (f"; because {nr['n_zero_replicates']} replicate(s) were zero, a small shift of {shift:.3g} was added before transforming" if shift > 0 else "")
             + ". The strain phenotype is the mean of the transformed replicates and the strain noise used for the wisam weights is their variance.\n")
    nfig = FIGS / "proteins" / f"{protein}_normalization.png"
    if nfig.exists():
        L.append(figure(nfig, f"{protein}: raw replicates, Box-Cox likelihood profile and transformed replicates.", [nfig, nfig.with_suffix(".pdf")]))
    # --- heritability + transcript pairing
    L.append("## Heritability\n")
    txt = f"Broad-sense heritability of the protein is **{f2(hr['H2'])}** (95% CI {f2(hr['CI_low'])}–{f2(hr['CI_high'])}; BH-adjusted P {fp(hr['p_value_BH'])}). "
    if pr["status"] == "no_pair":
        txt += f"No transcript pair is defined for this protein ({pr['note']})."
    else:
        shared = " This is a shared multi-gene probe, so the transcript signal is not gene-specific." if pr["status"] == "paired_shared_probe" else ""
        sym = short_symbol(pr["probe_symbol"])
        txt += (f"Its paired transcript is **{sym}** (array probe `{pr['Transcript_probe']}`, paired from the measured peptide `{pr['peptide']}`). "
                f"Transcript H² is **{f2(hr['tx_H2'])}** (95% CI {f2(hr['tx_low'])}–{f2(hr['tx_high'])}; BH-adjusted P {fp(hr['tx_p_BH'])}).{shared}")
    L.append(txt + "\n")
    # --- scan
    L.append("# pQTL genome scan\n")
    thr = num(ld.get("Threshold_0.95"))
    L.append(f"Genome scan of the Box-Cox phenotype with wisam weights (miqtl `scan.h2lmm`, 11 imputations); significance thresholds come from 1000 permutations "
             f"(GEV fit to the permutation maxima of −log10 P). Strongest locus: chr{ld.get('CH', '—')}:{f2(ld.get('Mb'))} Mb, P = {fp(ld.get('P.Val'))}, "
             f"adjusted P = {fp(ld.get('Adjusted_P.val'))}; 95% threshold = {f2(thr)} (−log10 P).\n")
    sp = OUT / "figs" / protein / "full" / "scans_and_adjusted0.95.pdf"
    if sp.exists():
        png = pdf_to_png(sp, PP / protein / f"{protein}_genome_scan.png", page=1)
        L.append(figure(png, f"{protein} genome-wide pQTL scan with the permutation threshold.", [stage(sp, protein, "genome_scan.pdf"), stage(OUT / "figs" / protein / "full" / "scans_0.95.pdf", protein, "scan_by_chromosome.pdf")]))
    # --- significant peaks
    L.append("# Significant peak(s)\n")
    if peaks:
        for pk in peaks:
            L.append(locus_blocks(protein, pk))
    else:
        L.append("No genome-wide significant pQTL for this protein (adjusted P ≥ 0.05), so credible intervals, Merge, TIMBR and mediation were not run.\n")
    # --- zero-inflated subset
    if protein in ZERO_INFLATED:
        L.append("# Zero-inflated subset scan\n")
        L.append(f"{protein} has strains with no detectable protein (all replicates zero). A second scan was run on the strains with measurable variation "
                 "(normalized within-strain variance > 0.01), re-normalized with its own Box-Cox λ, with its own permutation threshold. "
                 "The binary (detected vs. not detected) scan has not been re-run with the new normalization.\n")
        zp = ZOUT / "figs" / protein / "sub" / "scans_and_adjusted0.95.pdf"
        if zp.exists():
            png = pdf_to_png(zp, PP / protein / f"{protein}_sub_scan.png", page=1)
            L.append(figure(png, f"{protein} zero-inflated subset genome scan.", [stage(zp, protein, "sub_scan.pdf")]))
        zpk = sig_for(protein, zsig)
        if zpk:
            for pk in zpk:
                L.append(locus_blocks(protein, pk))
        else:
            L.append("The subset scan has no genome-wide significant peak.\n")
    (ROOT / f"{protein}.qmd").write_text("\n".join(L))


# --------------------------------------------------------------------------- integrated table / csv
def integrated_rows():
    rows = []
    for p in PROTEINS:
        k = p.upper()
        h, n, l, pr = her[k], norm[k], lead.get(k, {}), pairs[k]
        peaks = sig_for(p)
        types = sorted({cis.get((k, pk["Locus"]), {}).get("pqtl_type", "") for pk in peaks} - {""})
        rows.append(dict(
            Protein=k, lambda_=n["lambda"], lead_chr=l.get("CH", ""), lead_Mb=l.get("Mb", ""), lead_P=l.get("P.Val", ""), lead_adj_P=l.get("Adjusted_P.val", ""),
            n_significant_peaks=len(peaks), pqtl_type="/".join(types) if types else "",
            protein_H2=h["H2"], protein_H2_CI_low=h["CI_low"], protein_H2_CI_high=h["CI_high"], protein_H2_BH_P=h["p_value_BH"],
            transcript=(pr["probe_symbol"] if pr["status"] != "no_pair" else ""), transcript_pair_status=pr["status"],
            transcript_H2=h["tx_H2"], transcript_H2_CI_low=h["tx_low"], transcript_H2_CI_high=h["tx_high"], transcript_H2_BH_P=h["tx_p_BH"]))
    return rows


def overview_table(rows):
    t = ["| Protein | λ | Strongest pQTL (chr:Mb) | P | Adjusted P | Significant peaks (type) | Protein H² (95% CI) | Paired transcript | Transcript H² (95% CI) | Transcript BH P |",
         "|---|---:|---|---:|---:|---|---:|---|---:|---:|"]
    for r in rows:
        h2 = f"{f2(r['protein_H2'])} ({f2(r['protein_H2_CI_low'])}, {f2(r['protein_H2_CI_high'])})"
        th2 = (f"{f2(r['transcript_H2'])} ({f2(r['transcript_H2_CI_low'])}, {f2(r['transcript_H2_CI_high'])})" if num(r["transcript_H2"]) is not None else "—")
        shared = " †" if r["transcript_pair_status"] == "paired_shared_probe" else ""
        sigtxt = f"{r['n_significant_peaks']} ({r['pqtl_type']})" if r["n_significant_peaks"] else "0"
        t.append(f"| {r['Protein']} | {f2(r['lambda_'])} | {r['lead_chr']}:{f2(r['lead_Mb'])} | {fp(r['lead_P'])} | {fp(r['lead_adj_P'])} | {sigtxt} | {h2} | "
                 f"{md(short_symbol(r['transcript']) if r['transcript'] else '')}{shared} | {th2} | {fp(r['transcript_H2_BH_P'])} |")
    return "\n".join(t)


# --------------------------------------------------------------------------- pages
def followup_table():
    f = OUT / "pqtl_summary_table.csv"
    if not f.exists():
        return "_The locus follow-up summary (top Merge SDP, TIMBR allelic series, mediators) is added once Merge has finished._\n"
    shutil.copy2(f, FIGS / "pqtl_followup_summary.csv")
    t = ["### Locus follow-up summary\n", f"Per significant locus: the top Merge SDP (founder-haplotype split), the highest-posterior TIMBR allelic series at the peak, and the transcripts for which mediation is the top-supported model. {dl_line('Download (CSV)', FIGS / 'pqtl_followup_summary.csv')}.\n",
         "| Protein | Peak (chr:Mb) | pQTL type | Top SDP | Top TIMBR allelic series | Mediators |", "|---|---|---|---|---|---|"]
    for r in load_csv(f):
        c = cis.get((r["Protein"].upper(), r["Locus"]), {})
        t.append(f"| {r['Protein'].upper()} | {r['CH']}:{f2(r['Mb'])} | {c.get('pqtl_type') or '—'} | {md(r.get('top_sdp'))} | {md(r.get('top_timbr_allelic_series'))} | {md(r.get('mediators'))} |")
    return "\n".join(t) + "\n"


def write_overview(rows):
    csvp = FIGS / "integrated_summary.csv"
    with csvp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for src, dst in ((OUT / "pqtl_cis_trans.csv", "pqtl_cis_trans.csv"), (ZOUT / "pqtl_cis_trans.csv", "pqtl_cis_trans_zeroinflated.csv"),
                     (OUT / "protein_transcript_pairs.csv", "protein_transcript_pairs.csv")):
        if src.exists():
            shutil.copy2(src, FIGS / dst)
    for name in ("figure3_pqtl_scan_values.csv", "figure3_pqtl_significant_peaks.csv"):
        if (OUT / "figs" / "all" / name).exists():
            shutil.copy2(OUT / "figs" / "all" / name, FIGS / name)
    ov = sorted(p.stem.rsplit("chr", 1)[1] for p in FIGS.glob("figure3_pqtl_overlay_chr*.png"))
    n_sig = len(sig)
    sig_prot = sorted({r["Protein"] for r in sig})
    ct = sorted(cis.items())
    sig_rows = ["| Protein | Peak (chr:Mb) | Adjusted P | 90% CI (Mb) | Encoding gene | pQTL type |", "|---|---|---:|---|---|---|"]
    for pk in sig + zsig:
        c = cis.get((pk["Protein"].upper(), pk["Locus"]), {})
        sig_rows.append(f"| {pk['Protein'].upper()}{' (subset scan)' if pk.get('Subset') == 'sub' else ''} | {pk['CH']}:{pk['Mb']:.2f} | {fp(pk['Adj_P'])} | {pk['CI_start']:.1f}–{pk['CI_end']:.1f} | {md(c.get('pqtl_gene'))} | **{c.get('pqtl_type') or '—'}** |")
    eqtl = RNA_TAG
    page = f'''# Results overview {{.unnumbered}}

This is the cross-protein overview of the **{TAG}** analysis: protein phenotypes were Box-Cox normalized on the individual replicates, scanned for pQTL with wisam-weighted miqtl models, and followed up at each significant locus with credible intervals, Merge, TIMBR and RNA mediation. Protein symbols are uppercase; **CES2 is excluded** because of a protein QC issue. Every figure is shown here and can be downloaded with the link beneath it; per-protein results are in the protein chapters.

## Cross-protein summary

One row per protein: normalization λ, the strongest pQTL, significant peaks with their cis/trans type, protein broad-sense H², and the paired transcript with its H². A † marks a transcript measured by a shared multi-gene probe. {dl_line('Download the table (CSV)', csvp)}.

{overview_table(rows)}

## Significant pQTL

{len(sig)} genome-wide significant peaks in {len(sig_prot)} proteins (plus {len(zsig)} in the zero-inflated subset scan). A peak is **cis** when the gene encoding the protein lies within its 90% credible interval and **trans** otherwise ({dl_line('CSV', FIGS / 'pqtl_cis_trans.csv')}). Credible intervals, gene tracks and follow-up for each are in the [credible-interval appendix](credible_intervals.qmd) and the protein chapters.

{chr(10).join(sig_rows)}

## Heritability and strain structure

{tfig('figure2_heritability_profiles_split', 'Broad-sense heritability of each protein and its paired transcript (left) and clustered Box-Cox protein profiles across strains (right).', extra=('figure2_heritability_protein_and_transcript.csv',))}
Protein H² is the strain share of replicate variance (random-intercept REML on the Box-Cox replicates) with parametric-bootstrap 95% intervals; transcript H² uses individual-sample RINT values with scan date as a fixed effect. † marks a shared multi-gene probe and — a protein without a defensible transcript pair. Alternative layouts: {dl_line('single panel (PNG)', FIGS / 'figure2_heritability_profiles_blackgreen.png')} · {dl_line('combined colours (PNG)', FIGS / 'figure2_heritability_profiles.png')}.

### Protein–transcript pairing

{tfig('protein_transcript_correlation_heatmap', 'Correlation across CC strains between each protein and its paired transcript (outlined diagonal) and the other paired transcripts.', extra=('protein_transcript_correlation_matrix.csv',))}
Proteins were paired to transcripts from the measured peptide, not the protein name ({dl_line('pairing table (CSV)', FIGS / 'protein_transcript_pairs.csv')}). Proteins with no pair (CES1, CYP2C9) are omitted.

### Protein correlations

{tfig('protein_genetic_phenotypic_correlation_heatmap', 'Trait correlations: H² on the diagonal, genetic correlation (squares) below, phenotypic correlation (circles) above.', extra=('protein_genetic_correlation_lower_matrix.csv',))}
{tfig('protein_phenotype_correlation_heatmap', 'Phenotypic correlation (off-diagonal) with the genetic variance of each z-scored trait on the diagonal.', extra=('protein_phenotype_correlation_matrix.csv',))}
Genetic covariance was estimated on z-scored traits (method of moments on CC-line means, REML genetic variances), so genetic variances are on a common scale; correlations and H² are unaffected by the scaling. Matrices: {dl_line('genetic covariance (CSV)', FIGS / 'protein_genetic_correlation_H2_diagonal_genetic_covariance.csv')} · {dl_line('genetic correlation (CSV)', FIGS / 'protein_genetic_correlation_H2_diagonal_genetic_correlation.csv')}.

{tfig('transcript_genetic_correlation_H2_diagonal', 'Coding-transcript genetic correlations with broad-sense H² on the diagonal.')}

## Genome-wide QTL overview

### pQTL

{tfig('figure3_pqtl_genome_heatmap', 'Genome-wide pQTL scan signal (-log10 P) for every protein; circles mark significant peaks.', extra=('figure3_pqtl_scan_values.csv', 'figure3_pqtl_significant_peaks.csv'))}
Chromosomes with significant peaks for two or more proteins are shown as overlays (line opacity scales with the chromosome-wide signal) and as zoomed heatmaps.

{chr(10).join(tfig(f'figure3_pqtl_overlay_chr{c}', f'pQTL overlay, chromosome {c}.') for c in ov)}
{chr(10).join(tfig(f'figure3_pqtl_zoom_chr{c}', f'pQTL signal across proteins, chromosome {c}.') for c in ov)}

### eQTL

{tfig('figure3_eqtl_location_heatmap', 'eQTL peak position (x) against transcript gene position (y); the diagonal is cis.', extra=('figure3_eqtl_location_heatmap_data.csv',))}
The eQTL map is built from the corrected RNA data ({eqtl}); it does not depend on the protein normalization.

## Protein variation and transcriptome context

### Protein clustering

{tfig('protein_pca_biplot', 'PCA biplot of the protein profiles across strains, coloured by gap-statistic cluster.', extra=('protein_pca_loadings.csv', 'protein_pca_strain_scores.csv'))}
{tfig('protein_cluster_gapstat', 'Gap statistic for the number of strain clusters.', width='60%', extra=('protein_cluster_gapstat_values.csv',))}
{tfig('protein_cluster_heatmap', 'Strains clustered on their z-scored protein profiles.', extra=('protein_cluster_assignments.csv', 'protein_cluster_kw_drivers.csv'))}

### Transcriptome

{tfig('transcriptome_rint_hierarchical_heatmap', 'RINT-normalized transcriptome, hierarchically clustered.', extra=('transcriptome_rint_variance.csv',))}
{tfig('transcriptome_rint_all_transcript_binned_correlation', 'Binned correlation structure across the RINT-normalized transcriptome.', extra=('transcriptome_rint_binned_correlation.csv',))}
The correlation overview bins transcripts by expression variance for display; it is not the full 20,666 × 20,666 matrix.

## Follow-up coverage

Mediation, TIMBR and Merge were run at every significant locus; what is currently available is itemised in [Analysis progress](status.qmd).

{followup_table()}'''
    (ROOT / "overview.qmd").write_text(page)


def write_ci_appendix():
    L = ["# Credible intervals {.unnumbered}\n",
         "Every genome-wide significant pQTL locus with its 90% bootstrap credible interval (100 Bayesian bootstrap replicates of the chromosome scan; the 80% interval is boxed in the gene-track plots). "
         "A locus is **cis** when the gene encoding the protein lies inside the interval. Follow-up figures (Merge, TIMBR, mediation) are in each protein's chapter.\n",
         "| Protein | Chr | Peak (Mb) | Adjusted P | 90% CI (Mb) | Width (Mb) | Encoding gene | pQTL type |\n|---|---|---:|---:|---|---:|---|---|"]
    allp = sig + zsig
    for pk in allp:
        c = cis.get((pk["Protein"].upper(), pk["Locus"]), {})
        L.append(f"| {pk['Protein'].upper()}{' (subset)' if pk.get('Subset') == 'sub' else ''} | {pk['CH']} | {pk['Mb']:.2f} | {fp(pk['Adj_P'])} | "
                 f"{pk['CI_start']:.2f}–{pk['CI_end']:.2f} | {pk['CI_end'] - pk['CI_start']:.1f} | {md(c.get('pqtl_gene'))} | **{c.get('pqtl_type') or '—'}** |")
    L.append("")
    for pk in allp:
        protein = pk["Protein"]
        subset = pk.get("Subset", "full")
        fdir = fig_dir(subset) / protein / ("sub" if subset == "sub" else "full")
        L.append(f"## {protein.upper()}, chromosome {pk['CH']}" + (" (zero-inflated subset scan)" if subset == "sub" else "") + "\n")
        cg = fdir / f"ci_genes_chr{pk['CH']}.png"
        if cg.exists():
            L.append(figure(stage(cg, protein, f"ci_genes_{subset}_chr{pk['CH']}.png"), f"{protein.upper()}: genes within the chromosome {pk['CH']} credible interval.",
                            [stage(cg, protein, f"ci_genes_{subset}_chr{pk['CH']}.png")]))
        L.append(f"Interval {pk['CI_start']:.2f}–{pk['CI_end']:.2f} Mb around the peak at {pk['Mb']:.2f} Mb. Follow-up: see the [{protein.upper()} chapter]({protein}.qmd).\n")
    (ROOT / "credible_intervals.qmd").write_text("\n".join(L))


def write_status(rows):
    full_peaks = sig
    def state_counts(kind):
        have = 0
        for pk in full_peaks:
            c = pk["CH"]
            fd = OUT / "figs" / pk["Protein"] / "full"
            if kind == "merge":
                have += merge_state(pk["Protein"], pk["Locus"]) == "complete"
            elif kind == "mergeplot":
                have += (fd / f"mergeplot_{pk['Locus']}.png").exists()
            elif kind == "timbr_sweep":
                have += (OUT / pk["Protein"] / "TIMBR_ci" / f"ch{c}_res_list.RData").exists() or (OUT / pk["Protein"] / "TIMBR_ci" / f"ch{c}_singlelocus_res.RData").exists()
            elif kind == "timbr_plot":
                t = OUT / "figs" / pk["Protein"] / "TIMBR_ci"
                have += bool(list(t.glob(f"*chr{c}*"))) if t.exists() else 0
            elif kind == "mediation":
                have += (fd / f"targeted_mediation_chr{c}_posterior_bars.pdf").exists()
        return have
    n = len(full_peaks)
    scans = sum((OUT / "figs" / p / "full" / "scans_and_adjusted0.95.pdf").exists() for p in PROTEINS)
    t = [f"# Analysis progress {{.unnumbered}}\n",
         f"This page tracks the **{TAG}** analysis. CES2 is excluded because of a protein QC issue. Counts are computed from the files on disk when this page is built.\n",
         "| Analysis | Status | Note |", "|---|---|---|",
         f"| Replicate-level Box-Cox normalization | {len(PROTEINS)}/{len(PROTEINS)} proteins | λ per protein from `boxcox(y ~ strain)`; small shift only for proteins with zeros (Cyp2c39, Cyp2c50, Mrp3, Oct1.2). |",
         f"| pQTL genome scans (11 imputations, 1000 permutations) | {scans}/{len(PROTEINS)} proteins | wisam weights from the Box-Cox-scale replicate variance. |",
         f"| Significant pQTL peaks with 90% credible intervals | {n} peaks in {len({r['Protein'] for r in sig})} proteins | plus {len(zsig)} in the zero-inflated subset scan; all classified cis/trans. |",
         "| Zero-inflated subset scans (Cyp2c39, Cyp2c50) | 2/2 | the binary detected/not-detected scan has not been re-run with the new normalization. |",
         f"| Merge (diplotype fine-mapping) | {state_counts('merge')}/{n} full-set loci complete | {state_counts('mergeplot')}/{n} plots available; wide credible intervals take many hours. |",
         f"| TIMBR credible-interval sweep | {state_counts('timbr_sweep')}/{n} loci | {state_counts('timbr_plot')}/{n} loci have plots (plots need the Merge profile). |",
         f"| RNA mediation (bmediatR) | {state_counts('mediation')}/{n} loci | scan weights and phenotype come from the same tag. |",
         "| Protein broad-sense H² | 26/26 | 95% parametric-bootstrap intervals. |",
         f"| Transcript pairing and H² | {sum(r['transcript_pair_status'] != 'no_pair' for r in rows)}/{len(PROTEINS)} proteins paired | peptide-based pairing; CES1 and CYP2C9 have no defensible pair. |",
         "| Genetic covariance / correlation | Available | z-scored traits; point estimates (no uncertainty intervals). |",
         "| Protein clustering | Available | gap statistic, PCA, heatmap. |",
         "| eQTL map and transcriptome overview | Available | independent of the protein normalization. |",
         "",
         "## Pending",
         "- Merge, TIMBR plots and mediation overlays using the Merge-derived plotting window are completed once the running Merge jobs finish; rebuild the book afterwards (`python3 build_book_repbc.py && quarto render`).",
         "- The binary (detected vs not detected) scans for the zero-inflated proteins use the earlier definition."]
    (ROOT / "status.qmd").write_text("\n".join(t) + "\n")


def write_index_and_methods():
    (ROOT / "index.qmd").write_text(f'''# Preface {{.unnumbered}}

This book collects results from a pharmacogenomics screen of 26 drug-metabolizing enzyme and transporter (DMET) proteins measured across Collaborative Cross (CC) mouse strains (a 27th, CES2, is excluded for a protein QC issue). It reports the **{TAG}** analysis: protein abundances were Box-Cox normalized on the individual replicates, scanned for protein QTL (pQTL), and followed up at each significant locus.

The [results overview](overview.qmd) has the cross-protein summary, heritability, correlations, genome-wide QTL maps and clustering. [Analysis progress](status.qmd) lists what has been completed, and the [credible-interval appendix](credible_intervals.qmd) has every significant locus with its cis/trans call. The protein chapters give per-protein detail:

- **Normalization** -- raw, Box-Cox likelihood profile and transformed replicates.
- **Heritability** -- broad-sense H² of the protein and of its paired transcript.
- **pQTL genome scan** -- `miqtl` scan with permutation-derived thresholds.
- **Significant peaks** -- credible interval, cis/trans, Merge fine-mapping, TIMBR allelic series and RNA mediation.

Every figure is displayed in the page and has a download link (PNG, PDF or CSV) beneath it. See [Methods](intro.qmd) for how the data were generated and analysed.
''')
    (ROOT / "intro.qmd").write_text(f'''# Methods {{.unnumbered}}

## Proteomics

Raw protein-abundance data were generated using protein isolation, proteolytic digestion, and quantitative targeted absolute proteomics (QTAP) analysis, as previously described in Qasem et al. (2020). Peptides were selected to avoid confounding from genetic variation where possible. For **Cyp2c50**, **Cyp3a11** and **Oatp1a4** the peptide contains a missense variant in some Collaborative Cross founder/strain haplotypes. Three to four replicate measurements were made per protein per strain (45 CC strains).

## Normalization

For each protein a Box-Cox exponent λ was estimated on the individual replicates (`boxcox(y ~ strain)`, grid −2 to 2) with **no +1 shift**; a small shift (half the smallest positive replicate) is added only for proteins that contain zeros (Cyp2c39, Cyp2c50, Mrp3, Oct1.2). Every replicate is transformed with λ; the strain phenotype is the mean of the transformed replicates and the strain noise is their variance. For zero-inflated proteins (Cyp2c39, Cyp2c50) a second **subset** scan keeps strains whose normalized within-strain variance exceeds 0.01 and re-estimates λ on that subset.

## Genome scans

pQTL scans use [`miqtl`](https://cran.r-project.org/package=miqtl) (`scan.h2lmm`, haplotype-dosage model with 11 multiple imputations) on the strain phenotype, weighted by wisam precision weights computed from the transformed-scale noise. Genome-wide significance thresholds come from 1000 permutations (6 imputations each) with a generalized extreme value fit to the permutation maxima of −log10 P (95th percentile); adjusted P-values use the same fit. Credible intervals are 90% intervals from 100 Bayesian bootstrap replicates of the chromosome scan (80% intervals are drawn in the gene-track plots).

## cis and trans pQTL

A significant peak is **cis** if the gene encoding the protein lies inside its credible interval (same chromosome, mm10 gene body overlapping the interval) and **trans** otherwise. The encoding gene comes from the peptide-based pairing below (e.g. Ent1 → Slc29a1); CYP2C9 has no mouse ortholog and receives no call.

## Protein–transcript pairing

Each protein was paired to a transcript from the measured peptide (exact match against mouse UniProt sequences for the relevant gene families), not from the protein name. Where several array probes cover the gene, a single-gene probe is preferred, then the Affymetrix suffix (_at, _a_at, _s_at, _x_at), then highest expression. Proteins with no defensible pair (CES1: no Ces probes on the array; CYP2C9: no mouse ortholog) have no transcript; shared multi-gene probes (UGT1A1, UGT1A6, UGT2A2, CYP2C50) are marked † because their signal is not gene-specific.

## Heritability and genetic correlation

Broad-sense heritability is Vg/(Vg + Ve) from a strain random-intercept REML model on the Box-Cox replicates with a 200-replicate parametric-bootstrap interval; transcript H² uses individual-sample RINT values with scan date as a fixed effect. **Genetic covariance is computed on z-scored traits**: each protein's replicate-level values are z-scored before the REML fit and before the CC-line means used for off-diagonal covariances, so genetic variances share one scale; genetic correlations and H² are unaffected. Off-diagonal covariances use method of moments on CC-line means (protein replicates are not animal-matched across traits, so residual cross-trait covariance is assumed zero) and the correlation matrix is projected to the nearest positive-semidefinite matrix. Phenotypic correlations are Pearson correlations of the strain means.

## RNA mediation

Mediation analysis uses `bmediatR` (Bayesian mediation analysis for QTL) to test whether transcript abundance mediates a protein QTL, with the scan's own model weights, using corrected RNA data (duplicate probes resolved, all 20,666 transcripts).

## TIMBR

TIMBR (Tree-based Inference of Multiple-allele Regression) infers whether the founder haplotypes at a QTL collapse into fewer functional alleles, using a Chinese Restaurant Process prior over haplotype groupings with a phylogenetic backbone. It is run across the full credible interval with the scan's model weights.

## Merge / diplotype fine-mapping

At each significant locus, diplotype-resolved ("merge") analysis re-tests association at the level of imputed founder-haplotype dosage across the credible interval, with the scan's model weights, to narrow the candidate region and flag significant SNPs.

## Clustering

Strains are clustered (Ward.D2, Euclidean) on the z-scored Box-Cox protein profiles; the number of clusters is chosen with the gap statistic.

::: {{.callout-note}}
Full bibliographic details for Qasem et al. (2020) and Oreper et al. (2017) are not yet recorded in this project; please add the complete citations to `references.bib`.
:::
''')


def copy_book_assets():
    for name in ("transcriptome_rint_hierarchical_heatmap", "transcriptome_rint_all_transcript_binned_correlation"):
        for ext in ("png", "pdf"):
            src = ROOT / "figures" / "rerun_2026-09-24_v2" / f"{name}.{ext}"
            if src.exists() and not (FIGS / f"{name}.{ext}").exists():
                shutil.copy2(src, FIGS / f"{name}.{ext}")


def main():
    copy_book_assets()
    rows = integrated_rows()
    for p in PROTEINS:
        build_chapter(p)
    write_overview(rows)
    write_ci_appendix()
    write_status(rows)
    write_index_and_methods()
    print(f"Wrote index, intro (Methods), overview, status, credible_intervals and {len(PROTEINS)} protein chapters for {TAG}")


if __name__ == "__main__":
    main()
