"""Compatibility wrapper for running ARGWAS from a source checkout.

Prefer installing the package and running the ``argwas`` command, but this file
keeps the old workflow available from the repository root:

    python main.py -o out -p input.ped -t trees -a ancestral.txt
"""

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from argwas.cli import main


if __name__ == "__main__":
    main()
