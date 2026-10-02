#!/usr/bin/env python3
"""Point protein chapter genome-scan sections to the v2 scan PDFs."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCAN_DIR = ROOT / "figures" / "rerun_2026-09-24_v2" / "pqtl_per_protein"


def main():
    updated = 0
    missing = []
    for chapter in ROOT.glob("*.qmd"):
        if chapter.stem.lower() in {"ces2"}:
            continue
        text = chapter.read_text()
        start = text.find("## pQTL Genome Scans")
        end = text.find("## Peak Locus", start)
        if start < 0 or end < 0:
            continue
        pdf = next((p for p in SCAN_DIR.glob("*.pdf") if p.stem.lower() == chapter.stem.lower()), None)
        heading = "## pQTL Genome Scans (2026-09-24 v2; Box–Cox)"
        if pdf:
            relative = pdf.relative_to(ROOT).as_posix()
            section = (f"{heading}\n\n"
                       f"[Download the v2 genome scan PDF]({relative}).\n\n"
                       f"<iframe src=\"{relative}\" width=\"100%\" height=\"760px\" "
                       "title=\"2026-09-24 v2 pQTL genome scan\"></iframe>\n\n")
        else:
            missing.append(chapter.stem)
            section = (f"{heading}\n\n"
                       "No per-protein genome scan PDF was present in the saved 2026-09-24 v2 outputs. "
                       "This chapter does not display a scan from an older run in its place.\n\n")
        chapter.write_text(text[:start] + section + text[end:])
        updated += 1
    print(f"Updated v2 scan section in {updated} protein chapters")
    if missing:
        print("No saved v2 scan PDF for: " + ", ".join(missing))


if __name__ == "__main__":
    main()
