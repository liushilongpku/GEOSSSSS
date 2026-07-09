#!/usr/bin/env python3

import argparse
import base64
import json
import math
import os
import shlex
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_HOST = "server1"
DEFAULT_REMOTE_DIR = "/home/hello/codes/20260508/EXAMPLE/T10"
DEFAULT_PVDS = {
  "xita75": "reservoir_mefl_xita75.pvd",
  "xita80": "reservoir_mefl_xita80.pvd",
  "xita90": "reservoir_mefl_xita90.pvd",
}


REMOTE_ANALYSIS_SCRIPT = r"""
import json
import math
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


SECONDS_PER_YEAR = 31536000.0


def fracture_files_from_vtm(vtm_path):
    tree = ET.parse(vtm_path)
    root = tree.getroot()
    files = []
    for dataset in root.iter("DataSet"):
        rel_path = dataset.attrib.get("file", "")
        if "/Fracture/" in rel_path and rel_path.endswith(".vtu"):
            files.append(os.path.join(os.path.dirname(vtm_path), rel_path))
    return files


def analyze_dataset(base_dir, pvd_name, max_year):
    pvd_path = os.path.join(base_dir, pvd_name)
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    rows = []

    for dataset in root.iter("DataSet"):
        timestep_seconds = float(dataset.attrib["timestep"])
        year = timestep_seconds / SECONDS_PER_YEAR
        if year > max_year + 1e-9:
            continue

        vtm_path = os.path.join(base_dir, dataset.attrib["file"])
        projections = []
        count_state3 = 0

        for vtu_path in fracture_files_from_vtm(vtm_path):
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(vtu_path)
            reader.Update()
            grid = reader.GetOutput()
            cell_data = grid.GetCellData()

            fracture_state = cell_data.GetArray("fractureState")
            displacement_jump = cell_data.GetArray("displacementJump")
            normal_vector = cell_data.GetArray("normalVector")
            if fracture_state is None or displacement_jump is None or normal_vector is None:
                continue

            fracture_state_np = vtk_to_numpy(fracture_state)
            displacement_jump_np = vtk_to_numpy(displacement_jump)
            normal_vector_np = vtk_to_numpy(normal_vector)

            mask = fracture_state_np == 3
            if np.any(mask):
                normal_projection = np.einsum(
                    "ij,ij->i",
                    displacement_jump_np[mask],
                    normal_vector_np[mask],
                )
                projections.append(normal_projection)
                count_state3 += int(mask.sum())

        mean_normal_jump = math.nan
        if projections:
            mean_normal_jump = float(np.concatenate(projections).mean())

        rows.append(
            {
                "year": int(round(year)),
                "timestep_seconds": int(round(timestep_seconds)),
                "mean_normal_jump": mean_normal_jump,
                "count_state3": count_state3,
            }
        )

    return rows


def main():
    request = json.load(sys.stdin)
    base_dir = request["remote_dir"]
    max_year = float(request["max_year"])
    datasets = request["datasets"]

    output = {}
    for label, pvd_name in datasets.items():
        output[label] = analyze_dataset(base_dir, pvd_name, max_year)

    json.dump(output, sys.stdout, allow_nan=True)


if __name__ == "__main__":
    main()
"""


def parse_args():
  parser = argparse.ArgumentParser(
    description=(
      "Analyze remote GEOS fracture outputs and compare the yearly mean of "
      "displacementJump along the normal direction for fractureState == 3 cells."
    )
  )
  parser.add_argument("--host", default=DEFAULT_HOST, help="SSH host or alias.")
  parser.add_argument(
    "--remote-dir",
    default=DEFAULT_REMOTE_DIR,
    help="Remote directory containing the PVD/VTM/VTU files.",
  )
  parser.add_argument(
    "--max-year",
    type=float,
    default=30.0,
    help="Maximum year to include in the analysis.",
  )
  parser.add_argument(
    "--output-dir",
    default="analysis_outputs/remote_fracture_jump",
    help="Local directory for CSV and figure outputs.",
  )
  parser.add_argument(
    "--dataset",
    action="append",
    default=[],
    metavar="LABEL=PVD",
    help=(
      "Dataset mapping. May be passed multiple times. "
      "Defaults to xita75/xita80/xita90 if omitted."
    ),
  )
  return parser.parse_args()


def parse_dataset_args(dataset_args):
  if not dataset_args:
    return DEFAULT_PVDS.copy()

  datasets = {}
  for item in dataset_args:
    if "=" not in item:
      raise ValueError(f"Invalid --dataset value: {item!r}")
    label, pvd_name = item.split("=", 1)
    label = label.strip()
    pvd_name = pvd_name.strip()
    if not label or not pvd_name:
      raise ValueError(f"Invalid --dataset value: {item!r}")
    datasets[label] = pvd_name
  return datasets


def run_remote_analysis(host, remote_dir, datasets, max_year):
  request = {
    "remote_dir": remote_dir,
    "datasets": datasets,
    "max_year": max_year,
  }
  remote_input = json.dumps(request)
  encoded_script = base64.b64encode(REMOTE_ANALYSIS_SCRIPT.encode("utf-8")).decode("ascii")
  remote_command = (
    "import base64; "
    f"exec(base64.b64decode({encoded_script!r}).decode('utf-8'))"
  )
  command = ["ssh", host, f"python3 -c {shlex.quote(remote_command)}"]
  completed = subprocess.run(
    command,
    input=remote_input,
    text=True,
    capture_output=True,
    check=False,
  )
  if completed.returncode != 0:
    raise RuntimeError(
      "Remote analysis failed.\n"
      f"stdout:\n{completed.stdout}\n"
      f"stderr:\n{completed.stderr}"
    )

  stdout = completed.stdout.strip()
  if not stdout:
    raise RuntimeError("Remote analysis returned no data.")
  return json.loads(stdout)


def build_dataframe(results):
  rows = []
  for label, entries in results.items():
    for entry in entries:
      rows.append(
        {
          "dataset": label,
          "year": entry["year"],
          "timestep_seconds": entry["timestep_seconds"],
          "mean_normal_jump": entry["mean_normal_jump"],
          "count_state3": entry["count_state3"],
        }
      )

  frame = pd.DataFrame(rows)
  frame.sort_values(["dataset", "year"], inplace=True)
  return frame


def write_outputs(frame, output_dir):
  output_dir.mkdir(parents=True, exist_ok=True)

  csv_path = output_dir / "fracture_state3_normal_displacementjump_yearly.csv"
  frame.to_csv(csv_path, index=False)

  plt.style.use("seaborn-v0_8-whitegrid")
  fig, ax = plt.subplots(figsize=(10, 6))

  for dataset, group in frame.groupby("dataset"):
    group = group.sort_values("year")
    ax.plot(
      group["year"],
      group["mean_normal_jump"],
      marker="o",
      linewidth=2,
      markersize=4,
      label=dataset,
    )

  ax.set_xlabel("Year")
  ax.set_ylabel("Mean displacementJump · normalVector")
  ax.set_title("FractureState = 3 Cells: Yearly Mean Normal displacementJump")
  ax.legend()
  fig.tight_layout()

  plot_path = output_dir / "fracture_state3_normal_displacementjump_yearly.png"
  fig.savefig(plot_path, dpi=200)
  plt.close(fig)

  return csv_path, plot_path


def summarize(frame):
  lines = []
  for dataset, group in frame.groupby("dataset"):
    valid = group[group["count_state3"] > 0].sort_values("year")
    if valid.empty:
      lines.append(f"{dataset}: no fractureState == 3 cells within the selected years")
      continue

    first = valid.iloc[0]
    last = valid.iloc[-1]
    lines.append(
      (
        f"{dataset}: first active year = {int(first['year'])}, "
        f"year {int(last['year'])} mean = {last['mean_normal_jump']:.6e}, "
        f"state3 cells at year {int(last['year'])} = {int(last['count_state3'])}"
      )
    )
  return "\n".join(lines)


def main():
  args = parse_args()
  datasets = parse_dataset_args(args.dataset)
  results = run_remote_analysis(args.host, args.remote_dir, datasets, args.max_year)
  frame = build_dataframe(results)
  output_dir = Path(args.output_dir).resolve()
  csv_path, plot_path = write_outputs(frame, output_dir)

  print(summarize(frame))
  print(f"CSV: {csv_path}")
  print(f"Plot: {plot_path}")


if __name__ == "__main__":
  main()
