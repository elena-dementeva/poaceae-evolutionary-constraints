#!/usr/bin/env python3

import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import argparse

# --- Аргументы командной строки ---
parser = argparse.ArgumentParser(description="Построение phyloP-графиков по чанкам данных.")
parser.add_argument("input_file", help="Путь к входному файлу CSV/TSV без заголовков")
parser.add_argument("-o", "--output_dir", default="phyloP_plots", help="Папка для сохранения графиков")
parser.add_argument("--chunk_size", type=int, default=200_000, help="Размер чанка (строк читаемых из файла)")
parser.add_argument("--highlight", nargs="+", default=["gene"], help="Типы фичей для точек")
parser.add_argument("--background", nargs="+", default=["CDS"], help="Типы фичей для фона")
parser.add_argument("--point_size", type=int, default=10, help="Размер точек")
parser.add_argument("--alpha", type=float, default=0.6, help="Прозрачность точек")
parser.add_argument("--background_alpha", type=float, default=0.4, help="Прозрачность фона")
parser.add_argument("--delimiter", default="\t", help="Разделитель входного файла")
parser.add_argument("--columns", nargs="+", default=[
    "chrom", "start", "end", "position", "phyloP", "p_value", "q_value",
    "significant", "constraint_type", "gff_chrom", "gff_start", "gff_end", "feature_type", "gff_attributes"
], help="Названия колонок (по умолчанию merged_with_all.tsv)")

args = parser.parse_args()

# --- Подготовка ---
os.makedirs(args.output_dir, exist_ok=True)

cmap_points = colormaps.get_cmap('tab10')
cmap_background = colormaps.get_cmap('Dark2')
default_color = 'lightgray'

# --- Чтение и обработка чанков ---
reader = pd.read_csv(args.input_file, sep=args.delimiter, header=None, chunksize=args.chunk_size)

for chunk_idx, chunk in enumerate(reader, 1):
    print(f"\n[INFO] Обработка чанка {chunk_idx}...")
    # Назначение колонок
    chunk.columns = args.columns
    chunk['feature_type'] = chunk['feature_type'].fillna('.')

    # Определяем окно по чанкy
    min_pos = int(chunk['end'].min())
    max_pos = int(chunk['end'].max())
    print(f"[INFO] Диапазон позиций: {min_pos} - {max_pos}")

    highlight_types_unique = list(set(chunk['feature_type']) & set(args.highlight))
    background_types_unique = list(set(chunk['feature_type']) & set(args.background))

    color_map = {
        ftype: cmap_points(i / len(highlight_types_unique))
        for i, ftype in enumerate(highlight_types_unique)
    }
    background_color_map = {
        ftype: cmap_background(i / len(background_types_unique))
        for i, ftype in enumerate(background_types_unique)
    }

    # Группировка фичей по позиции
    types_per_pos = chunk.groupby("end")["feature_type"].apply(set).to_dict()
    highlight_positions = {
        pos for pos, types in types_per_pos.items()
        if types & set(args.highlight)
    }

    print("Построение графика для чанка")
    # --- Построение графика для чанка ---
    fig, ax = plt.subplots(figsize=(14, 4))

    # Полная ось X
    x_full = np.arange(min_pos, max_pos + 1)
    ax.plot(x_full, [np.nan] * len(x_full), color='none')

    # Подготовка данных
    unique_positions_df = chunk.drop_duplicates(subset=['end'])
    pos_array = unique_positions_df['end'].values
    phyloP_array = unique_positions_df['phyloP'].values

    colors = [
        color_map[sorted(types_per_pos.get(pos) & set(args.highlight))[0]]
        if pos in highlight_positions and types_per_pos.get(pos) & set(args.highlight)
        else default_color
        for pos in pos_array
    ]

        # Рисуем точки
    for x, y, c in zip(pos_array, phyloP_array, colors):
        ax.vlines(x, ymin=0, ymax=y, color=c, linewidth=args.point_size * 0.1, alpha=args.alpha)

    # Фоновые области
    background_df = chunk[chunk['feature_type'].isin(args.background)].drop_duplicates(
        subset=["gff_start", "gff_end", "feature_type"]
    )
    for _, row in background_df.iterrows():
        ftype = row['feature_type']
        b_color = background_color_map.get(ftype, 'orange')
        ax.axvspan(row["gff_start"], row["gff_end"], color=b_color, alpha=args.background_alpha)

    # Оформление осей
    ax.set_xlim(min_pos, max_pos)
    ax.set_ylim(chunk["phyloP"].min() - 0.5, chunk["phyloP"].max() + 0.5)
    ax.set_xlabel("Genomic position")
    ax.set_ylabel("phyloP -log10(p-value)")
    ax.set_title(f"PhyloP: positions {min_pos}-{max_pos}")
    ax.grid(True, linestyle="--", alpha=0.3)

    # Легенда
    legend_elements = []
    for ftype in highlight_types_unique:
        legend_elements.append(
            plt.Line2D([0], [0], color=color_map[ftype], lw=4, label=f'{ftype} (point)')
        )
    for ftype in background_types_unique:
        legend_elements.append(
            plt.Line2D([0], [0], color=background_color_map[ftype], lw=10, alpha=args.background_alpha, label=f"{ftype} (bg)")
        )
    if legend_elements:
        ax.legend(handles=legend_elements, title="Feature Types", bbox_to_anchor=(1.01, 1), loc='upper left')

    plt.tight_layout()
    plot_path = os.path.join(args.output_dir, f"phyloP_chunk_{chunk_idx}.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Сохранён график: {plot_path}")
