#!/usr/bin/env python3
"""Run the SPE6 7x7x8 fine-grid restart chain locally.

Key information: outputs default to /tmp, completed segments are reused, and
each segment advances exactly one 0.25-year SPE6 reporting interval.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPOSITORY = next(parent for parent in ROOT.parents if (parent / ".git").exists())
DEFAULT_GEOS = REPOSITORY / "build/bin/geosx"
DEFAULT_OUTPUT = Path("/tmp/fine_grid_spe6_7x7x8_saturated_tol5e2_max12")


def manifest() -> list[dict[str, str]]:
    """Read generated segment records while skipping traceability comments."""

    lines = [
        line
        for line in (ROOT / "segment_manifest.csv").read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]
    return list(csv.DictReader(lines))


def final_restart(output_dir: Path, stem: str) -> Path | None:
    """Return the latest restart root without its .root suffix."""

    roots = sorted(output_dir.glob(f"{stem}_restart_*.root"))
    return roots[-1].with_suffix("") if roots else None


def completed(output_dir: Path, stem: str) -> bool:
    """Require both GEOS' finish marker and a restart file."""

    log = output_dir / "run.log"
    return bool(
        final_restart(output_dir, stem)
        and log.is_file()
        and "Finished at" in log.read_text(encoding="utf-8", errors="replace")
    )


def main() -> None:
    """Run or reuse the requested prefix of the 20-segment chain."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--geos", type=Path, default=DEFAULT_GEOS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--segments", type=int, default=20, help="number of 0.25-year segments")
    args = parser.parse_args()
    geos = args.geos.resolve()
    output_root = args.output.resolve()
    if not geos.is_file():
        raise FileNotFoundError(geos)
    records = manifest()
    if not 1 <= args.segments <= len(records):
        raise ValueError(f"segments must be in [1, {len(records)}]")

    previous_restart: Path | None = None
    for row in records[: args.segments]:
        deck = ROOT / row["deck"]
        stem = deck.stem
        output_dir = output_root / stem
        if completed(output_dir, stem):
            previous_restart = final_restart(output_dir, stem)
            print(f"reuse {stem}", flush=True)
            continue
        if output_dir.exists() and any(output_dir.iterdir()):
            attempt = 1
            while (archive := output_dir.with_name(f"{stem}_failed_attempt_{attempt:02d}")).exists():
                attempt += 1
            output_dir.rename(archive)
            print(f"preserve incomplete output: {archive}", flush=True)
        output_dir.mkdir(parents=True, exist_ok=True)
        command = [str(geos), "-i", str(deck), "-o", str(output_dir)]
        if previous_restart is not None:
            command.extend(("-r", str(previous_restart)))
        print(f"run {stem}", flush=True)
        with (output_dir / "run.log").open("w", encoding="utf-8") as log:
            result = subprocess.run(
                command,
                cwd=REPOSITORY,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if result.returncode != 0 or not completed(output_dir, stem):
            raise RuntimeError(f"GEOS failed for {stem}; inspect {output_dir / 'run.log'}")
        previous_restart = final_restart(output_dir, stem)


if __name__ == "__main__":
    main()
