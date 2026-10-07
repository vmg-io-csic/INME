#!/usr/bin/env python3
"""Build the figures in one command, on any platform.

    python make_figures.py                 draw what has already been exported
    python make_figures.py --run           run every experiment first
    python make_figures.py ReversePhi      draw one figure

This is a convenience wrapper. The two steps it runs work on their own:

    matlab -batch "run_all_experiments"
    python plot_figures.py

Running the experiments needs MATLAB. Drawing does not: if the export folder
is in the repository, a clone can reproduce every figure with Python alone.
"""

import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main(argv):
    args = list(argv)
    if "--run" in args:
        args.remove("--run")
        matlab = shutil.which("matlab")
        if matlab is None:
            sys.exit("matlab is not on PATH; load it first, or drop --run to "
                     "draw from the export that is already there.")
        r = subprocess.run([matlab, "-batch", "run_all_experiments"], cwd=HERE)
        if r.returncode != 0:
            sys.exit(f"the experiments exited with code {r.returncode}")

    r = subprocess.run([sys.executable, str(HERE / "plot_figures.py")] + args,
                       cwd=HERE)
    sys.exit(r.returncode)


if __name__ == "__main__":
    main(sys.argv[1:])
