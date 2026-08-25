#!/usr/bin/env python3
"""Plot oil-component mass recovery for the four P0 control-variable runs."""

# Purpose: compare mass-based oil-component recovery without mixing it with the SPE6 volume-recovery reference.
# Key information: component 0 in compAmount is used, and each curve is normalized by its own initial mass.

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


YEAR = 31_536_000.0


def load_mass_recovery(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        time = handle["compAmount Time"][:].ravel() / YEAR
        oil_component_mass = handle["compAmount"][:][:, 0, 0]
    recovery = (oil_component_mass[0] - oil_component_mass) / oil_component_mass[0] * 100.0
    return time, recovery


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old", type=Path, required=True)
    parser.add_argument("--oldcomp-scaledpc", type=Path, required=True)
    parser.add_argument("--livecomp-unscaledpc", type=Path, required=True)
    parser.add_argument("--revised", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    runs = [
        ("旧版：无初始溶解气，Pc 未缩放", args.old, "#0072B2", "o"),
        ("控制 A：无初始溶解气，Pc 缩放", args.oldcomp_scaledpc, "#009E73", "D"),
        ("控制 B：活油初始化，Pc 未缩放", args.livecomp_unscaledpc, "#CC79A7", "x"),
        ("修正版：活油初始化，Pc 缩放", args.revised, "#D55E00", "^"),
    ]

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
    for label, path, color, marker in runs:
        time, recovery = load_mass_recovery(path)
        ax.plot(time, recovery, color=color, linewidth=2.2, marker=marker, markersize=4.2,
                label=f"{label}（{recovery[-1]:.2f}%）")

    ax.set_xlim(0.0, 5.05)
    ax.set_ylim(0.0, 66.0)
    ax.set_xlabel("时间 / 年")
    ax.set_ylabel("基质油组分质量采收率 / %")
    ax.set_title("P0 SPE6 单块：油组分质量采收率")
    ax.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
    ax.legend(loc="lower right", frameon=True, framealpha=0.95)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
