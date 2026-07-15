#!/usr/bin/env python3
"""
Plot arbitrary GEOS DPDP Mandel pressure histories against the analytical
solution.

Examples:
  python plot_GEOS_vs_analytical.py \
    --case "FIM 0.911=/path/to/fim_0911" \
    --case "FIM 1.0=/path/to/fim_1000" \
    --case "Seq=/path/to/seq_full" \
    --out /path/to/pressure_comparison.png

  python plot_GEOS_vs_analytical.py \
    --search-root inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs \
    --out inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/pressure_comparison.png

Each case directory must contain:
  pressure_matrix_history.hdf5
  pressure_fracture_history.hdf5
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "geos_mandel_mplconfig"))

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
VALIDATION_DIR = SCRIPT_DIR.parent
DEFAULT_RUNS_DIR = VALIDATION_DIR / "runs"

sys.path.insert(0, str(SCRIPT_DIR))
import dpdp_mandel_analytical as an

PM0 = 4.55e5
PF0 = 4.88e5
REQUIRED_HISTORY_FILES = ("pressure_matrix_history.hdf5", "pressure_fracture_history.hdf5")


@dataclass(frozen=True)
class Case:
  label: str
  output_dir: Path


def existing_path(value: str) -> Path:
  path = Path(value).expanduser().resolve()
  if not path.exists():
    raise argparse.ArgumentTypeError(f"path does not exist: {path}")
  return path


def has_pressure_histories(path: Path) -> bool:
  return path.is_dir() and all((path / name).is_file() for name in REQUIRED_HISTORY_FILES)


def center_index(element_centers: np.ndarray) -> int:
  """Element nearest the x-min no-flow axis and mid-height z=0.015."""
  x = element_centers[:, 0]
  z = element_centers[:, 2]
  return int(np.argmin(np.abs(x - x.min()) + np.abs(z - 0.015)))


def load_pressure(output_dir: Path, filename: str) -> tuple[np.ndarray, np.ndarray]:
  path = output_dir / filename
  with h5py.File(path, "r") as h5:
    time = np.asarray(h5["pressure Time"])[:, 0]
    pressure = np.asarray(h5["pressure"])
    element_centers = np.asarray(h5["pressure elementCenter"])[0]
  index = center_index(element_centers)
  mask = time > 0.0
  return time[mask] / an.t0, pressure[mask, index]


def read_xml_text(output_dir: Path) -> str:
  chunks: list[str] = []
  for path in sorted(output_dir.glob("*.xml")):
    try:
      chunks.append(path.read_text(errors="replace"))
    except OSError:
      pass
  return "\n".join(chunks)


def compact_label(text: str) -> str:
  text = re.sub(r"[_-]+", " ", text)
  text = re.sub(r"\s+", " ", text).strip()
  return text or "GEOS case"


def infer_label(output_dir: Path) -> str:
  """Infer a readable label from copied XML metadata or the directory name."""
  xml_text = read_xml_text(output_dir)
  lower = f"{output_dir.name}\n{xml_text}".lower()

  if "couplingtype=\"sequential\"" in lower or re.search(r"\bseq\b|sequential", lower):
    return f"Seq: {compact_label(output_dir.name)}"

  scale = re.search(r"crossstorageoffdiagscale=\"([^\"]+)\"", lower)
  if "couplingtype=\"fullyimplicit\"" in lower or re.search(r"\bfim\b|fullyimplicit", lower):
    if scale:
      return f"FIM scale={scale.group(1)}: {compact_label(output_dir.name)}"
    return f"FIM: {compact_label(output_dir.name)}"

  return compact_label(output_dir.name)


def parse_case(value: str) -> Case:
  if "=" in value:
    label, raw_path = value.split("=", 1)
    label = label.strip()
    if not label:
      raise argparse.ArgumentTypeError("--case label cannot be empty")
  else:
    raw_path = value
    label = ""

  output_dir = Path(raw_path).expanduser().resolve()
  if not output_dir.exists():
    raise argparse.ArgumentTypeError(f"case output directory does not exist: {output_dir}")
  if not has_pressure_histories(output_dir):
    raise argparse.ArgumentTypeError(
        f"case output directory does not contain required pressure histories: {output_dir}"
    )
  return Case(label or infer_label(output_dir), output_dir)


def discover_cases(search_roots: list[Path]) -> list[Case]:
  candidates: list[Path] = []
  for root in search_roots:
    if has_pressure_histories(root):
      candidates.append(root)
    if root.exists():
      candidates.extend(path for path in root.rglob("*") if has_pressure_histories(path))

  unique: dict[Path, Case] = {}
  for path in sorted(set(candidates), key=lambda item: str(item)):
    unique[path] = Case(infer_label(path), path)
  return list(unique.values())


def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(
      description="Plot any number of DPDP Mandel GEOS pressure histories against the analytical solution."
  )
  parser.add_argument(
      "--case",
      action="append",
      type=parse_case,
      default=[],
      metavar="LABEL=OUTPUT_DIR",
      help="Add one GEOS output case. May be repeated. LABEL= is optional.",
  )
  parser.add_argument(
      "--search-root",
      action="append",
      type=existing_path,
      default=[],
      help="Search recursively for GEOS output directories. May be repeated.",
  )
  parser.add_argument(
      "--include-default-runs",
      action="store_true",
      help="Also search ../runs relative to this script, when it exists.",
  )
  parser.add_argument("--out", type=Path, help="Output PNG path. A PDF and summary CSV are written beside it.")
  parser.add_argument("--tau-min", type=float, default=1e-4, help="Minimum dimensionless time shown on the plot.")
  parser.add_argument("--tau-max", type=float, default=3e4, help="Maximum dimensionless time shown on the plot.")
  parser.add_argument(
      "--analytical-points",
      type=int,
      default=150,
      help="Number of analytical inversion points. Larger values are smoother but slower.",
  )
  return parser.parse_args()


def deduplicate_cases(cases: list[Case]) -> list[Case]:
  by_path: dict[Path, Case] = {}
  for case in cases:
    if case.output_dir not in by_path:
      by_path[case.output_dir] = case
  return list(by_path.values())


def log_interp(x: np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
  mask = (x > 0.0) & (x >= xp.min()) & (x <= xp.max())
  out = np.full_like(x, np.nan, dtype=float)
  out[mask] = np.interp(np.log10(x[mask]), np.log10(xp), fp)
  return out


def choose_output_path(args: argparse.Namespace, cases: list[Case]) -> Path:
  if args.out:
    return args.out.expanduser().resolve()

  if cases:
    common = Path(os.path.commonpath([str(case.output_dir) for case in cases]))
    if common.is_file():
      common = common.parent
    return common / "geos_vs_analytical_pressure.png"

  return Path.cwd() / "geos_vs_analytical_pressure.png"


def style_for(index: int) -> tuple[object, str, str]:
  color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
  color = color_cycle[index % len(color_cycle)] if color_cycle else None
  linestyles = ["-", "--", "-.", ":"]
  markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]
  return color, linestyles[index % len(linestyles)], markers[index % len(markers)]


def plot_pressure_cases(cases: list[Case], args: argparse.Namespace) -> tuple[Path, Path, Path]:
  sol = an.solve(tau_min=min(1e-5, args.tau_min), tau_max=args.tau_max, n=args.analytical_points)
  analytical = {
      "matrix": (sol["tau"], sol["pm"]),
      "fracture": (sol["tau"], sol["pf"]),
  }

  output_path = choose_output_path(args, cases)
  output_path.parent.mkdir(parents=True, exist_ok=True)
  pdf_path = output_path.with_suffix(".pdf")
  csv_path = output_path.with_name(output_path.stem + "_summary.csv")

  fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), sharex=True)
  axes[0].semilogx(sol["tau"], sol["pm"], color="black", lw=2.0, label="Analytical")
  axes[1].semilogx(sol["tau"], sol["pf"], color="black", lw=2.0, label="Analytical")

  summary_rows: list[dict[str, object]] = []
  legend_handles = [Line2D([0], [0], color="black", lw=2.0, label="Analytical")]

  for index, case in enumerate(cases):
    color, linestyle, marker = style_for(index)
    tau_m, pm = load_pressure(case.output_dir, "pressure_matrix_history.hdf5")
    tau_f, pf = load_pressure(case.output_dir, "pressure_fracture_history.hdf5")

    legend_handles.append(Line2D([0], [0], color=color, lw=1.55, ls=linestyle, marker=marker,
                                 ms=3.0, label=case.label))

    for ax, medium, tau, pressure, p0 in (
        (axes[0], "matrix", tau_m, pm, PM0),
        (axes[1], "fracture", tau_f, pf, PF0),
    ):
      normalized = pressure / p0
      marker_every = max(1, len(tau) // 45)
      ax.semilogx(
          tau,
          normalized,
          color=color,
          ls=linestyle,
          lw=1.55,
          marker=marker,
          ms=2.2,
          markevery=marker_every,
      )

      analytical_tau, analytical_y = analytical[medium]
      analytical_at_geos = log_interp(tau, analytical_tau, analytical_y)
      valid = ~np.isnan(analytical_at_geos)
      peak_index = int(np.argmax(normalized))
      summary_rows.append({
          "case": case.label,
          "medium": medium,
          "output_dir": str(case.output_dir),
          "n_samples": len(normalized),
          "peak_tau": tau[peak_index],
          "peak_p_over_p0": normalized[peak_index],
          "final_tau": tau[-1],
          "final_p_over_p0": normalized[-1],
          "mean_abs_error_vs_analytical": float(np.mean(np.abs(normalized[valid] - analytical_at_geos[valid]))),
          "max_abs_error_vs_analytical": float(np.max(np.abs(normalized[valid] - analytical_at_geos[valid]))),
      })

  for ax, title in zip(
      axes,
      ("Matrix pressure at center/no-flow axis", "Fracture pressure at center/no-flow axis"),
  ):
    ax.set_title(title)
    ax.set_xlabel(r"Dimensionless time $\tau=t/t_0$" + f"  ($t_0={an.t0:.3g}$ s)")
    ax.set_ylabel(r"$p/p_0^+$")
    ax.set_xlim(args.tau_min, args.tau_max)
    ax.set_ylim(0.0, 1.16)
    ax.grid(True, which="both", alpha=0.28)

  legend_columns = 1 if len(legend_handles) <= 5 else 2
  axes[1].legend(handles=legend_handles, frameon=False, fontsize=8, loc="best", ncol=legend_columns)

  fig.suptitle("DPDP Mandel-Cryer pressure validation: GEOS cases vs analytical", fontsize=12)
  fig.tight_layout(rect=(0, 0, 1, 0.94))
  fig.savefig(output_path, dpi=180)
  fig.savefig(pdf_path)
  plt.close(fig)

  with csv_path.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0].keys()))
    writer.writeheader()
    writer.writerows(summary_rows)

  return output_path, pdf_path, csv_path


def main() -> int:
  args = parse_args()
  if args.tau_min <= 0 or args.tau_max <= args.tau_min:
    raise ValueError("--tau-min must be positive and --tau-max must be larger than --tau-min")
  if args.analytical_points < 8:
    raise ValueError("--analytical-points must be at least 8")

  search_roots = list(args.search_root)
  if args.include_default_runs or not args.case:
    if DEFAULT_RUNS_DIR.exists():
      search_roots.append(DEFAULT_RUNS_DIR)

  cases = deduplicate_cases(list(args.case) + discover_cases(search_roots))
  if not cases:
    print("No GEOS output cases were provided or discovered.", file=sys.stderr)
    print("Use --case LABEL=OUTPUT_DIR or --search-root RUNS_DIR.", file=sys.stderr)
    return 2

  print(f"t0 = {an.t0:.6g} s")
  for case in cases:
    print(f"{case.label}: {case.output_dir}")

  output_path, pdf_path, csv_path = plot_pressure_cases(cases, args)
  print(f"saved {output_path}")
  print(f"saved {pdf_path}")
  print(f"saved {csv_path}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
