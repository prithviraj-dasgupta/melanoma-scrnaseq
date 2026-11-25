#!/usr/bin/env bash
# Run the whole analysis from download to final figures (about 3-5 minutes on a laptop).
set -euo pipefail
cd "$(dirname "$0")/.."

uv run python scripts/01_download.py
uv run python scripts/02_inspect_data.py
uv run python scripts/03_qc.py
uv run python scripts/04_cluster.py
uv run python scripts/05_annotate.py
uv run python scripts/06_cnv.py
uv run python scripts/07_programs.py
