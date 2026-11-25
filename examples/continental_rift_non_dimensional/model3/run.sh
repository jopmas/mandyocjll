#!/usr/bin/env bash
# Creates the inputs (if needed) and runs the example: bash run.sh [number of MPI ranks]
set -e
cd "$(dirname "$0")"
PYTHON=${PYTHON:-python}
MANDYOCJLL=${MANDYOCJLL:-../../../bin/mandyocjll}
[ -f param.txt ] || "$PYTHON" make_inputs.py
"$MANDYOCJLL" -n "${1:-1}" -seed 0,2 -strain_seed 0.0,1.0
