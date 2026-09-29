#!/usr/bin/env python3
"""Backward-compatible entry point for the current v2 report builder."""

from pathlib import Path
import runpy


if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).with_name("build_v2_report.py")),
                   run_name="__main__")
