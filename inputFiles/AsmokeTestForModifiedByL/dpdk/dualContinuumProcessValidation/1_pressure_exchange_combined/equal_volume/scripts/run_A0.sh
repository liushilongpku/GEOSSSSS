#!/usr/bin/env bash
# Purpose: run the A0 pressure-exchange input into its case-local archive.
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 GEOSX [extra geosx args...]" >&2
  exit 2
fi
geosx="$1"
shift
case_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$case_root/runs/20260821"
exec "$geosx" -i "$case_root/A0_pressure_exchange.xml" -o "$case_root/runs/20260821" "$@"
