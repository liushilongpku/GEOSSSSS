#!/usr/bin/env bash
# Purpose: run the J0 local derivative check and archive its text result in this case.
set -euo pipefail

case_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$case_root/runs/20260821"
python3 "$case_root/scripts/check_J0.py" | tee "$case_root/runs/20260821/check_J0.txt"
