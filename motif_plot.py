#!/usr/bin/env python
"""
Cluster motif visualisation for GraphClust results.
Reads RESULTS/*/cluster.all files and produces a motif overview plot (PDF).

Usage:
  python motif_plot.py <results_dir> [output_file]
"""

import glob
import itertools
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import numpy as np
from matplotlib import pyplot as plt

# Colourblind-friendly categorical palette (Okabe-Ito inspired)
PALETTE = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
]


def parse_clusters(results_dir):
    """Read all cluster.all files from results_dir/*/cluster.all."""
    cluster_files = sorted(glob.glob(os.path.join(results_dir, "*", "cluster.all")))
    if not cluster_files:
        print("WARNING: no cluster.all files found in %s" % results_dir)
        return {}, {}, {}, {}

    palette = itertools.cycle(PALETTE)
    ranges = defaultdict(list)
    colors = defaultdict(list)
    orig_names = {}
    cluster_nums = {}

    for cluster_file in cluster_files:
        cluster_color = next(palette)
        with open(cluster_file) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) < 11 or parts[1] != "RESULT":
                    continue
                seq, start, end, strand = parts[0].split("#")
                ranges[seq].append((int(start), max(1, int(end) - int(start) + 1)))
                colors[seq].append(cluster_color)
                cluster_label = "cluster-%s" % parts[2]
                cluster_nums[cluster_label] = cluster_color
                if "ORIGHEAD" in parts:
                    idx = parts.index("ORIGHEAD")
                    if idx + 1 < len(parts):
                        orig_names[seq] = parts[idx + 1]
                elif seq not in orig_names:
                    orig_names[seq] = seq

    return ranges, colors, orig_names, cluster_nums


def plot_bar(ranges, colors, orig_names, cluster_nums, output_file):
    """Create horizontal bar chart of motif positions."""
    if not ranges:
        print("No data to plot.")
        return

    sorted_keys = sorted(ranges.keys())
    fig, ax = plt.subplots(figsize=(10, max(3, 0.4 * len(sorted_keys))))

    for i, k in enumerate(sorted_keys):
        ax.broken_barh(ranges[k], (i - 0.25, 0.5), facecolors=colors[k])

    ax.set_xlim(0)
    ax.set_xlabel("Position in sequence")
    ax.set_yticks(range(len(sorted_keys)))
    short_labels = []
    for k in sorted_keys:
        name = orig_names.get(k, k)
        if len(name) > 30:
            name = name[:27] + "..."
        short_labels.append("%s - %s" % (k, name))
    ax.set_yticklabels(short_labels, fontsize=7)
    ax.grid(True, alpha=0.3)

    patches = [
        mpatches.Patch(color=cluster_nums[lab], label=lab)
        for lab in sorted(cluster_nums)
    ]
    ax.legend(handles=patches, loc="upper right", fontsize=7, framealpha=0.8)

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()
    print("Motif plot saved to %s" % output_file)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python motif_plot.py <results_dir> [output_file]")
        sys.exit(1)

    results_dir = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else os.path.join(results_dir, "motif_plot.pdf")

    r, c, n, cn = parse_clusters(results_dir)
    plot_bar(r, c, n, cn, output_file)
