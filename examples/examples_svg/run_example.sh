#!/usr/bin/env bash
# Creates the inputs of an example from its SVG (tapioca-svg) and runs it:
#     bash run_example.sh NAME [number of MPI ranks]
# The files go to NAME/inputs. It needs tapIOca (conda environment tapioca); other options of tapioca-svg may be
# given in TAPIOCA_SVG_OPTIONS, e.g. TAPIOCA_SVG_OPTIONS="--dx 10 --dz 10".
set -e
cd "$(dirname "$0")"
export MANDYOCJLL=${MANDYOCJLL:-$(cd ../.. && pwd)/bin/mandyocjll}
[ -n "$1" ] && [ -d "$1" ] || { echo "Examples:"; find . -name "*.svg" ! -path "*/inputs/*" | sed -E 's|^\./(.*)/[^/]*$|    \1|' | sort; exit 1; }
cd "$1"
SVG=$(find . -maxdepth 1 -name "*.svg" | head -1)
[ -n "$SVG" ] || { echo "no SVG in $1"; exit 1; }
if [ ! -f inputs/param.txt ]; then
    tapioca-svg "$SVG" -o inputs --ncores "${2:-1}" --no-post-processing $TAPIOCA_SVG_OPTIONS
fi
cd inputs
bash run_*.sh
