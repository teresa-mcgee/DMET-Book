#!/usr/bin/env python3
"""
Convert PDFs to JPGs from a given input folder to an output folder.
Checks for blank PDFs and removes them.
If JPGs already exist, they are replaced by new ones.

Usage:
    python convert_pdfs.py <input_folder> <output_folder>
"""

import sys
from pathlib import Path
import fitz  # PyMuPDF
from PIL import Image, ImageChops

def is_blank_pdf(pdf_path: Path) -> bool:
    """Check if a PDF is blank (no text and no graphics)."""
    try:
        doc = fitz.open(pdf_path)
        if len(doc) == 0:
            return True
        for page in doc:
            if page.get_text().strip():
                return False
            pix = page.get_pixmap(dpi=30)
            img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
            bg = Image.new("RGB", img.size, (255, 255, 255))
            diff = ImageChops.difference(img, bg)
            if diff.getbbox() is not None:
                return False
        return True
    except Exception:
        return True

def pdf_to_jpgs(pdf_path: Path, output_folder: Path):
    """Convert all pages of a PDF into JPGs, replacing old files if they exist."""
    doc = fitz.open(pdf_path)
    pdf_name = pdf_path.stem
    subfolder = output_folder / pdf_name
    subfolder.mkdir(parents=True, exist_ok=True)

    for i, page in enumerate(doc, start=1):
        pix = page.get_pixmap(dpi=300)
        img_path = subfolder / f"{pdf_name}_page_{i:03d}.jpg"

        # If an old JPG exists, remove it before saving
        if img_path.exists():
            print(f"    🔄 Replacing existing {img_path.name}")
            img_path.unlink()

        pix.save(str(img_path))

    return subfolder

def main(input_folder: Path, output_folder: Path):
    if not input_folder.exists():
        print(f"❌ Input folder does not exist: {input_folder}")
        sys.exit(1)

    output_folder.mkdir(parents=True, exist_ok=True)
    print(f"Processing PDFs in: {input_folder}")
    print(f"Saving JPGs to: {output_folder}")

    for pdf_file in input_folder.glob("*.pdf"):
        if is_blank_pdf(pdf_file):
            print(f"  ⚠️ Removing blank PDF: {pdf_file.name}")
            pdf_file.unlink()
            continue

        print(f"  📄 Converting {pdf_file.name} ...")
        saved_to = pdf_to_jpgs(pdf_file, output_folder)
        print(f"  ✅ Saved images in {saved_to}")

    print("🎉 Conversion complete!")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_pdfs.py <input_folder> <output_folder>")
        sys.exit(1)

    input_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    main(input_dir, output_dir)
