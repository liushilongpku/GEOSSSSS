#!/usr/bin/env python3
"""Plot P0 recovery curves against the digitized SPE6 reference."""

# Purpose: compare old and revised P0 matrix oil recovery with the digitized SPE6 zero-fracture-Pc curve.
# Key information: recovery is computed from matrix oil saturation, with no smoothing or interpolation of simulation data.

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


YEAR = 31_536_000.0


def load_recovery(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        time = handle["phaseVolumeFraction Time"][:].ravel() / YEAR
        phase = handle["phaseVolumeFraction"][:][:, 0, :]
    recovery = (phase[0, 0] - phase[:, 0]) / phase[0, 0] * 100.0
    return time, recovery


def load_reference(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, delimiter=",", skiprows=3)
    return data[:, 0], data[:, 1] * 100.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old", type=Path, required=True, help="old P0 phase history")
    parser.add_argument("--revised", type=Path, required=True, help="revised P0 phase history")
    parser.add_argument("--oldcomp-scaledpc", type=Path, help="old composition with scaled Pc phase history")
    parser.add_argument("--livecomp-unscaledpc", type=Path, help="live composition with unscaled Pc phase history")
    parser.add_argument("--reference", type=Path, required=True, help="digitized SPE6 reference CSV")
    parser.add_argument("--output", type=Path, required=True, help="PNG output path")
    args = parser.parse_args()

    old_time, old_recovery = load_recovery(args.old)
    revised_time, revised_recovery = load_recovery(args.revised)
    reference_time, reference_recovery = load_reference(args.reference)
    oldcomp_scaled = load_recovery(args.oldcomp_scaledpc) if args.oldcomp_scaledpc else None
    livecomp_unscaled = load_recovery(args.livecomp_unscaledpc) if args.livecomp_unscaledpc else None

    cjk_font = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    if cjk_font.exists():
        font_manager.fontManager.addfont(cjk_font)
        cjk_family = font_manager.FontProperties(fname=cjk_font).get_name()
    else:
        cjk_family = "DejaVu Sans"
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [cjk_family, "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 11,
    })
    fig, ax = plt.subplots(figsize=(8.6, 5.6), layout="constrained")
    ax.plot(reference_time, reference_recovery, color="#222222", linestyle=(0, (4, 2)), linewidth=2.0,
            marker="s", markersize=4.5, label="SPE6数字化参考")
    ax.plot(old_time, old_recovery, color="#0072B2", linewidth=2.2, marker="o", markersize=4.2,
            label=f"旧版：无初始溶解气，Pc 未缩放（{old_recovery[-1]:.2f}%）")
    if oldcomp_scaled is not None:
        ax.plot(oldcomp_scaled[0], oldcomp_scaled[1], color="#009E73", linewidth=2.0, marker="D", markersize=3.8,
                label=f"控制 A：无初始溶解气，Pc 缩放（{oldcomp_scaled[1][-1]:.2f}%）")
    ax.plot(revised_time, revised_recovery, color="#D55E00", linewidth=2.2, marker="^", markersize=4.2,
            label=f"修正版：活油初始化，Pc 缩放（{revised_recovery[-1]:.2f}%）")
    if livecomp_unscaled is not None:
        ax.plot(livecomp_unscaled[0], livecomp_unscaled[1], color="#CC79A7", linewidth=2.0, marker="x", markersize=4.0,
                label=f"控制 B：活油初始化，Pc 未缩放（{livecomp_unscaled[1][-1]:.2f}%）")
    ax.set_xlim(0.0, 5.05)
    ax.set_ylim(0.0, 48.0)
    ax.set_xlabel("时间 / 年")
    ax.set_ylabel("基质油相体积采收率 / %")
    ax.set_title("P0 SPE6 单块：初始组分与毛管压力缩放影响")
    ax.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
    ax.legend(loc="lower right", frameon=True, framealpha=0.95)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
