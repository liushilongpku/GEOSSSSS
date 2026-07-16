#!/usr/bin/env python3
"""
Plot arbitrary GEOS DPDP Mandel pressure histories against the retained
Fig. 5c analytical pressure points.

Examples:
  python plot_GEOS_vs_analytical.py \
    --case "FIM 0.911=/path/to/fim_0911" \
    --case "FIM 1.0=/path/to/fim_1000" \
    --case "Seq=/path/to/seq_full" \
    --out /path/to/pressure_comparison.png \
    --show-script-analytical

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
DEFAULT_ANALYTICAL_DIR = VALIDATION_DIR / "analitical_result"
DEFAULT_MATRIX_REFERENCE = DEFAULT_ANALYTICAL_DIR / "fig5c_primary_analitical.csv"
DEFAULT_FRACTURE_REFERENCE = DEFAULT_ANALYTICAL_DIR / "fig5c_secondary_analitical.csv"

sys.path.insert(0, str(SCRIPT_DIR))
import dpdp_mandel_analytical as an

PM0 = 4.55e5
PF0 = 4.88e5
REQUIRED_HISTORY_FILES = ("pressure_matrix_history.hdf5", "pressure_fracture_history.hdf5")
DEFAULT_SAMPLE_TAU = np.asarray([1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1, 3, 10, 30, 100, 300, 1000, 3000])


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


def load_reference_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
  tau: list[float] = []
  value: list[float] = []
  with path.open(newline="") as stream:
    reader = csv.DictReader(stream)
    if not reader.fieldnames or "x" not in reader.fieldnames or "y" not in reader.fieldnames:
      raise ValueError(f"reference CSV must contain x,y columns: {path}")
    for row in reader:
      tau.append(float(row["x"]))
      value.append(float(row["y"]))
  if len(tau) < 2:
    raise ValueError(f"reference CSV must contain at least two rows: {path}")
  order = np.argsort(np.asarray(tau))
  return np.asarray(tau)[order], np.asarray(value)[order]


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
      description="Plot any number of DPDP Mandel GEOS pressure histories against the Fig. 5c CSV reference."
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
  parser.add_argument("--out", type=Path, help="Output PNG path. A PDF, summary CSV, and samples CSV are written beside it.")
  parser.add_argument("--tau-min", type=float, default=1e-4, help="Minimum dimensionless time shown on the plot.")
  parser.add_argument("--tau-max", type=float, default=3e4, help="Maximum dimensionless time shown on the plot.")
  parser.add_argument(
      "--matrix-reference",
      type=existing_path,
      default=DEFAULT_MATRIX_REFERENCE,
      help="Manual matrix/primary reference CSV with x,y columns.",
  )
  parser.add_argument(
      "--fracture-reference",
      type=existing_path,
      default=DEFAULT_FRACTURE_REFERENCE,
      help="Manual fracture/secondary reference CSV with x,y columns.",
  )
  parser.add_argument(
      "--show-script-analytical",
      action="store_true",
      help="Overlay the local analytical-script curve for diagnosis. Errors are still computed against the CSV reference.",
  )
  parser.add_argument(
      "--analytical-points",
      type=int,
      default=150,
      help="Number of analytical-script inversion points used only with --show-script-analytical.",
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


def finite_error_metrics(value: np.ndarray, reference: np.ndarray) -> tuple[float, float, int]:
  valid = np.isfinite(value) & np.isfinite(reference)
  if not np.any(valid):
    return float("nan"), float("nan"), 0
  error = np.abs(value[valid] - reference[valid])
  return float(np.mean(error)), float(np.max(error)), int(np.count_nonzero(valid))


def choose_output_path(args: argparse.Namespace, cases: list[Case]) -> Path:
  if args.out:
    return args.out.expanduser().resolve()

  if cases:
    common = Path(os.path.commonpath([str(case.output_dir) for case in cases]))
    if common.is_file():
      common = common.parent
    return common / "geos_vs_fig5c_pressure.png"

  return Path.cwd() / "geos_vs_fig5c_pressure.png"


def style_for(index: int) -> tuple[object, str, str]:
  color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
  color = color_cycle[index % len(color_cycle)] if color_cycle else None
  linestyles = ["-", "--", "-.", ":"]
  markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]
  return color, linestyles[index % len(linestyles)], markers[index % len(markers)]


def log_spaced_marker_indices(tau: np.ndarray, target_count: int = 45) -> list[int]:
  """Pick marker locations with roughly uniform spacing on the log-time axis."""
  valid = np.flatnonzero(np.isfinite(tau) & (tau > 0.0))
  if len(valid) <= target_count:
    return valid.tolist()

  log_tau = np.log10(tau[valid])
  targets = np.linspace(log_tau[0], log_tau[-1], target_count)
  picked: list[int] = []
  for target in targets:
    right = int(np.searchsorted(log_tau, target, side="left"))
    if right <= 0:
      chosen = 0
    elif right >= len(log_tau):
      chosen = len(log_tau) - 1
    else:
      left = right - 1
      chosen = left if abs(log_tau[left] - target) <= abs(log_tau[right] - target) else right
    index = int(valid[chosen])
    if not picked or picked[-1] != index:
      picked.append(index)
  return picked


def plot_pressure_cases(cases: list[Case], args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
  reference = {
      "matrix": load_reference_csv(args.matrix_reference),
      "fracture": load_reference_csv(args.fracture_reference),
  }
  script_solution = None
  if args.show_script_analytical:
    script_solution = an.solve(tau_min=min(1e-5, args.tau_min), tau_max=args.tau_max, n=args.analytical_points)

  output_path = choose_output_path(args, cases)
  output_path.parent.mkdir(parents=True, exist_ok=True)
  pdf_path = output_path.with_suffix(".pdf")
  csv_path = output_path.with_name(output_path.stem + "_summary.csv")
  samples_path = output_path.with_name(output_path.stem + "_samples.csv")

  fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), sharex=True)
  axes[0].semilogx(reference["matrix"][0], reference["matrix"][1], color="black", lw=2.0, label="Fig. 5c CSV")
  axes[1].semilogx(reference["fracture"][0], reference["fracture"][1], color="black", lw=2.0, label="Fig. 5c CSV")
  if script_solution is not None:
    axes[0].semilogx(script_solution["tau"], script_solution["pm"], color="0.45", lw=1.4, ls=":",
                     label="Analytical script")
    axes[1].semilogx(script_solution["tau"], script_solution["pf"], color="0.45", lw=1.4, ls=":",
                     label="Analytical script")

  summary_rows: list[dict[str, object]] = []
  legend_handles = [Line2D([0], [0], color="black", lw=2.0, label="Fig. 5c CSV")]
  if script_solution is not None:
    legend_handles.append(Line2D([0], [0], color="0.45", lw=1.4, ls=":", label="Analytical script"))

  sample_rows: list[dict[str, object]] = []
  for tau in DEFAULT_SAMPLE_TAU:
    row: dict[str, object] = {"tau": tau}
    for medium in ("matrix", "fracture"):
      row[f"manual_{medium}"] = log_interp(np.asarray([tau]), reference[medium][0], reference[medium][1])[0]
    sample_rows.append(row)

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
      ax.semilogx(
          tau,
          normalized,
          color=color,
          ls=linestyle,
          lw=1.55,
          marker=marker,
          ms=2.2,
          markevery=log_spaced_marker_indices(tau),
      )

      reference_tau, reference_y = reference[medium]
      reference_at_geos = log_interp(tau, reference_tau, reference_y)
      mean_error, max_error, n_error_samples = finite_error_metrics(normalized, reference_at_geos)
      peak_index = int(np.argmax(normalized))
      summary_rows.append({
          "case": case.label,
          "medium": medium,
          "output_dir": str(case.output_dir),
          "n_samples": len(normalized),
          "n_error_samples": n_error_samples,
          "reference": "manual_fig5c_csv",
          "peak_tau": tau[peak_index],
          "peak_p_over_p0": normalized[peak_index],
          "final_tau": tau[-1],
          "final_p_over_p0": normalized[-1],
          "mean_abs_error_vs_manual_csv": mean_error,
          "max_abs_error_vs_manual_csv": max_error,
      })

      sample_at_geos = log_interp(DEFAULT_SAMPLE_TAU, tau, normalized)
      sample_reference = log_interp(DEFAULT_SAMPLE_TAU, reference_tau, reference_y)
      safe_label = re.sub(r"[^0-9A-Za-z]+", "_", case.label).strip("_").lower() or f"case_{index}"
      for row_index, row in enumerate(sample_rows):
        row[f"{safe_label}_{medium}"] = sample_at_geos[row_index]
        if np.isfinite(sample_at_geos[row_index]) and np.isfinite(sample_reference[row_index]):
          row[f"{safe_label}_{medium}_abs_error"] = abs(sample_at_geos[row_index] - sample_reference[row_index])
        else:
          row[f"{safe_label}_{medium}_abs_error"] = float("nan")

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

  fig.suptitle("DPDP Mandel-Cryer pressure validation: GEOS cases vs Fig. 5c CSV", fontsize=12)
  fig.tight_layout(rect=(0, 0, 1, 0.94))
  fig.savefig(output_path, dpi=180)
  fig.savefig(pdf_path)
  plt.close(fig)

  with csv_path.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0].keys()))
    writer.writeheader()
    writer.writerows(summary_rows)

  with samples_path.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(sample_rows[0].keys()))
    writer.writeheader()
    writer.writerows(sample_rows)

  return output_path, pdf_path, csv_path, samples_path


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

  print(f"matrix reference = {args.matrix_reference}")
  print(f"fracture reference = {args.fracture_reference}")
  output_path, pdf_path, csv_path, samples_path = plot_pressure_cases(cases, args)
  print(f"saved {output_path}")
  print(f"saved {pdf_path}")
  print(f"saved {csv_path}")
  print(f"saved {samples_path}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
