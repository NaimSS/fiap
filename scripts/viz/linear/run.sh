#!/bin/bash
# Runs visualize.py in its own venv (scripts/viz/linear/.venv), installing the
# pinned requirements on first use.
#
#   ./run.sh <instance.fi|.fi.diameter> <m> [options]    (./run.sh --help)

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

VENV="$SCRIPT_DIR/.venv"
PY="$VENV/bin/python"
if [[ ! -x "$PY" ]]; then
    echo "Creating venv in $VENV ..." >&2
    python3 -m venv "$VENV" || { echo "error: could not create $VENV" >&2; exit 1; }
fi
[[ -x "$PY" ]] || { echo "error: $PY missing after venv creation" >&2; exit 1; }
STAMP="$VENV/.requirements.stamp"
if [[ ! -f "$STAMP" || "$SCRIPT_DIR/requirements.txt" -nt "$STAMP" ]]; then
    echo "Installing pinned dependencies into $VENV ..." >&2
    "$PY" -m pip install --quiet -r "$SCRIPT_DIR/requirements.txt" \
        || { echo "error: pip install into $VENV failed" >&2; exit 1; }
    touch "$STAMP"
fi

"$PY" "$SCRIPT_DIR/visualize.py" "$@"
