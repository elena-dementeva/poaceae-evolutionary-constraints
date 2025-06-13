import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps

parser = argparse.ArgumentParser(description="Build phyloP plots by data chunks.")
parser.add_argument(
    "input_file",
    help="Path to input CSV/TSV file without headers",
)
parser.add_argument(
    "-o",
    "--output_dir",
    default="phyloP_plots",
    help="Directory to save plots",
)
parser.add_argument(
    "--chunk_size",
    type=int,
    default=200_000,
    help="Chunk size (rows per read)",
)
parser.add_argument(
    "--highlight",
    nargs="+",
    default=["gene"],
    help="Feature types for points",
)
parser.add_argument(
    "--background",
    nargs="+",
    default=["CDS"],
    help="Feature types for background",
)
parser.add_argument(
    "--point_size",
    type=int,
    default=10,
    help="Point size",
)
parser.add_argument(
    "--alpha",
    type=float,
    default=0.6,
    help="Point transparency",
)
parser.add_argument(
    "--background_alpha",
    type=float,
    default=0.4,
    help="Background transparency",
)
parser.add_argument(
    "--delimiter",
    default="\t",
    help="Input file delimiter",
)
parser.add_argument(
    "--columns",
    nargs="+",
    default=[
        "chrom",
        "start",
        "end",
        "position",
        "phyloP",
        "p_value",
        "q_value",
        "significant",
        "constraint_type",
        "gff_chrom",
        "gff_start",
        "gff_end",
        "feature_type",
        "gff_attributes",
    ],
    help="Column names (default for merged_with_all.tsv)",
)
args = parser.parse_args()


os.makedirs(args.output_dir, exist_ok=True)

CMAP_POINTS = colormaps.get_cmap("tab10")
CMAP_BACKGROUND = colormaps.get_cmap("Dark2")
DEFAULT_COLOR = "lightgray"


reader = pd.read_csv(
    args.input_file,
    sep=args.delimiter,
    header=None,
    chunksize=args.chunk_size,
)

for chunk_idx, chunk in enumerate(reader, start=1):
    print(f"[INFO] Processing chunk {chunk_idx}...")
    chunk.columns = args.columns
    chunk["feature_type"] = chunk["feature_type"].fillna(".")

    min_pos = int(chunk["end"].min())
    max_pos = int(chunk["end"].max())
    print(f"[INFO] Position range: {min_pos} - {max_pos}")

    highlight_types = list(set(chunk["feature_type"]) & set(args.highlight))
    background_types = list(set(chunk["feature_type"]) & set(args.background))

    color_map = {
        ftype: CMAP_POINTS(i / len(highlight_types))
        for i, ftype in enumerate(highlight_types)
    }
    background_color_map = {
        ftype: CMAP_BACKGROUND(i / len(background_types))
        for i, ftype in enumerate(background_types)
    }

    types_per_pos = chunk.groupby("end")["feature_type"].apply(set).to_dict()
    highlight_positions = {
        pos for pos, types in types_per_pos.items() if types & set(args.highlight)
    }

    print("[INFO] Building plot for chunk")
    fig, ax = plt.subplots(figsize=(14, 4))

    x_full = np.arange(min_pos, max_pos + 1)
    ax.plot(x_full, [np.nan] * len(x_full), color="none")

    unique_df = chunk.drop_duplicates(subset=["end"])
    pos_array = unique_df["end"].values
    phyloP_array = unique_df["phyloP"].values

    colors = [
        (
            color_map[sorted(types_per_pos[pos] & set(args.highlight))[0]]
            if pos in highlight_positions
            else DEFAULT_COLOR
        )
        for pos in pos_array
    ]

    for x, y, c in zip(pos_array, phyloP_array, colors):
        ax.vlines(
            x,
            ymin=0,
            ymax=y,
            color=c,
            linewidth=args.point_size * 0.1,
            alpha=args.alpha,
        )

    background_df = chunk[chunk["feature_type"].isin(args.background)].drop_duplicates(
        subset=["gff_start", "gff_end", "feature_type"]
    )
    for _, row in background_df.iterrows():
        ftype = row["feature_type"]
        b_color = background_color_map.get(ftype, "orange")
        ax.axvspan(
            row["gff_start"],
            row["gff_end"],
            color=b_color,
            alpha=args.background_alpha,
        )

    ax.set_xlim(min_pos, max_pos)
    ax.set_ylim(
        chunk["phyloP"].min() - 0.5,
        chunk["phyloP"].max() + 0.5,
    )
    ax.set_xlabel("Genomic position")
    ax.set_ylabel("phyloP -log10(p-value)")
    ax.set_title(f"PhyloP: positions {min_pos}-{max_pos}")
    ax.grid(True, linestyle="--", alpha=0.3)

    legend_elements = []
    for ftype in highlight_types:
        legend_elements.append(
            plt.Line2D(
                [0],
                [0],
                color=color_map[ftype],
                lw=4,
                label=f"{ftype} (point)",
            )
        )
    for ftype in background_types:
        legend_elements.append(
            plt.Line2D(
                [0],
                [0],
                color=background_color_map[ftype],
                lw=10,
                alpha=args.background_alpha,
                label=f"{ftype} (bg)",
            )
        )
    if legend_elements:
        ax.legend(
            handles=legend_elements,
            title="Feature Types",
            bbox_to_anchor=(1.01, 1),
            loc="upper left",
        )

    plt.tight_layout()
    plot_path = os.path.join(args.output_dir, f"phyloP_chunk_{chunk_idx}.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Saved plot: {plot_path}")
