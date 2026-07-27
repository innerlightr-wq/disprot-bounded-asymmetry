#!/usr/bin/env python3
"""Execute the complete analysis pipeline in dependency order.

    python run_all.py

Regenerates every table, figure and diagnostic used by the manuscript from the
four DisProt JSON exports in data/raw. Run order matters: 02 consumes the table
written by 01, and 03/06 consume what 02 writes.

Optional environment overrides (both honoured by every script):

    DISPROT_DATA   directory holding the four DisProt JSON exports
                   (default: <repo>/data/raw)
    DISPROT_OUT    directory to write results into
                   (default: <repo>/results)

Nothing in results/deposited_originals/ is ever written by this pipeline; those
files are archival copies from the published Zenodo record and are left
untouched.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"

SCRIPTS = [
    "01_build_tables.py",
    "02_asymmetry.py",
    "03_figures.py",
    "04_retrieval_log.py",
    "05_dphi_test.py",
    "06_fig4_column.py",
]


def main() -> int:
    missing = [s for s in SCRIPTS if not (SRC / s).exists()]
    if missing:
        print(f"ERROR: missing script(s) in {SRC}: {', '.join(missing)}",
              file=sys.stderr)
        return 1

    for script in SCRIPTS:
        print(f"Running {script}...", flush=True)
        proc = subprocess.run([sys.executable, str(SRC / script)],
                              cwd=ROOT)
        if proc.returncode != 0:
            print(f"\nERROR: {script} exited with code {proc.returncode}. "
                  f"Pipeline halted; later stages depend on this output.",
                  file=sys.stderr)
            return proc.returncode

    print("Pipeline completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
