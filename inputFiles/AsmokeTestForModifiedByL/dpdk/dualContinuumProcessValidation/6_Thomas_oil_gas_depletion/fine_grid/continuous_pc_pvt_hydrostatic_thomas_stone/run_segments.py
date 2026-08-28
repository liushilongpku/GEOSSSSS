#!/usr/bin/env python3
"""Run the Thomas Eq. (25) deck-approximation restart chain sequentially.

Key information: reuse requires both a GEOS finish marker and final restart;
incomplete attempts are preserved rather than overwritten.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[5]
DEFAULT_GEOS = REPOSITORY / "build/bin/geosx"


def read_manifest(variant: str) -> list[dict[str, str]]:
    """Load generated segment records for one approximation variant."""

    lines = [
        line
        for line in (ROOT / "segment_manifest.csv").read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]
    return [row for row in csv.DictReader(lines) if row["variant"] == variant]


def final_restart(output_dir: Path, stem: str) -> Path | None:
    """Return the final restart root path without the .root suffix."""

    roots = sorted(output_dir.glob(f"{stem}_restart_*.root"))
    return roots[-1].with_suffix("") if roots else None


def completed(output_dir: Path, stem: str) -> bool:
    """Check both GEOS completion and restart-file availability."""

    log = output_dir / "run.log"
    restart = final_restart(output_dir, stem)
    return bool(restart and log.exists() and "Finished at" in log.read_text(encoding="utf-8", errors="replace"))


def archive_incomplete(output_dir: Path) -> Path | None:
    """Preserve a failed output directory before rerunning cleanly."""

    if not output_dir.exists() or not any(output_dir.iterdir()):
        return None
    attempt = 1
    while True:
        archive = output_dir.with_name(f"{output_dir.name}_failed_attempt_{attempt:02d}")
        if not archive.exists():
            output_dir.rename(archive)
            return archive
        attempt += 1


def main() -> None:
    """Run the generated Thomas Eq. (25) right-endpoint chain."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", choices=("left", "right", "both"), default="both")
    parser.add_argument("--geos", type=Path, default=DEFAULT_GEOS)
    args = parser.parse_args()
    variants = ("left", "right") if args.variant == "both" else (args.variant,)
    geos = args.geos.resolve()
    if not geos.is_file():
        raise FileNotFoundError(geos)

    for variant in variants:
        previous_restart: Path | None = None
        for row in read_manifest(variant):
            deck = ROOT / row["deck"]
            stem = deck.stem
            output_dir = ROOT / "runs" / variant / stem
            if completed(output_dir, stem):
                previous_restart = final_restart(output_dir, stem)
                print(f"[{variant}] reuse {stem}: {previous_restart}", flush=True)
                continue

            archive = archive_incomplete(output_dir)
            if archive is not None:
                print(f"[{variant}] preserve incomplete output: {archive}", flush=True)
            output_dir.mkdir(parents=True, exist_ok=True)

            command = [str(geos), "-i", str(deck), "-o", str(output_dir)]
            if previous_restart is not None:
                command.extend(("-r", str(previous_restart)))
            print(f"[{variant}] run {stem}", flush=True)
            with (output_dir / "run.log").open("w", encoding="utf-8") as log:
                result = subprocess.run(command, cwd=REPOSITORY, stdout=log, stderr=subprocess.STDOUT, check=False)
            if result.returncode != 0 or not completed(output_dir, stem):
                raise RuntimeError(f"GEOS failed for {variant}/{stem}; inspect {output_dir / 'run.log'}")
            previous_restart = final_restart(output_dir, stem)


if __name__ == "__main__":
    main()
