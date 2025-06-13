import argparse
import logging
import multiprocessing as mp
import os
import re
import time
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
from statsmodels.stats.multitest import multipletests

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_header(line):
    """Parse a fixedStep header line and extract chromosome, start position, and step size."""
    chrom_match = re.search(r"chrom=(\S+)", line)
    pos_match = re.search(r"start=(\d+)", line)
    step_match = re.search(r"step=(\d+)", line)

    if chrom_match and pos_match:
        chrom = chrom_match.group(1)
        start = int(pos_match.group(1))
        step = int(step_match.group(1)) if step_match else 1
        return chrom, start, step
    return None, None, None


def analyze_file_structure(input_file):
    """
    Analyze the WIG file structure by identifying all fixedStep blocks.
    Returns a list of block metadata plus total line counts.
    """
    blocks = []
    total_lines = 0
    data_lines = 0

    with open(input_file, "r") as f:
        file_pos = 0
        current_block = None

        for line_num, line in enumerate(f):
            total_lines += 1

            if line.startswith("fixedStep"):
                data_lines += 1
                chrom, start, step = parse_header(line)

                if current_block:
                    current_block["end_pos"] = file_pos
                    current_block["line_count"] = line_num - current_block["start_line"]
                    blocks.append(current_block)

                current_block = {
                    "header": line.strip(),
                    "chrom": chrom,
                    "start": start,
                    "step": step,
                    "file_pos": file_pos,
                    "start_line": line_num,
                }
            else:
                # Non-header lines are data lines when inside a block
                if current_block is not None:
                    data_lines += 1

            file_pos += len(line)

        if current_block:
            current_block["end_pos"] = file_pos
            current_block["line_count"] = (
                line_num - current_block["start_line"] if "line_num" in locals() else 0
            )
            blocks.append(current_block)

    logger.info(
        f"File analysis: {len(blocks)} blocks, {total_lines} total lines, {data_lines} data lines"
    )
    return blocks, total_lines, data_lines


def read_block_data(input_file, block_info):
    """
    Read data for a single fixedStep block from the WIG file.

    Args:
        input_file: Path to the WIG file.
        block_info: Dictionary with block metadata from analyze_file_structure.

    Returns:
        Tuple (chrom, positions, scores) or None on error.
    """
    chrom = block_info["chrom"]
    start_pos = block_info["start"]
    step = block_info["step"]
    file_pos = block_info["file_pos"]

    try:
        with open(input_file, "r") as f:
            f.seek(file_pos)
            _ = f.readline()  # skip header line

            positions = []
            scores = []
            current_pos = start_pos

            while True:
                line = f.readline()
                if not line or line.startswith("fixedStep"):
                    break
                line = line.strip()
                if not line:
                    continue
                try:
                    score = float(line)
                    positions.append(current_pos)
                    scores.append(score)
                    current_pos += step
                except ValueError:
                    continue

            return chrom, positions, scores
    except Exception as e:
        logger.error(f"Error reading block at position {file_pos}: {e}")
        return None


def process_blocks(input_file, blocks_to_process, process_id=0):
    """
    Process a set of WIG file blocks in one worker.

    Args:
        input_file: Path to the WIG file.
        blocks_to_process: List of block metadata dictionaries.
        process_id: Worker ID for logging.

    Returns:
        List of result dictionaries for each position.
    """
    results = []

    for block in blocks_to_process:
        block_data = read_block_data(input_file, block)
        if not block_data:
            continue

        chrom, positions, scores = block_data
        positions_array = np.array(positions)
        scores_array = np.array(scores)

        for pos, scr in zip(positions_array, scores_array):
            status = (
                "constrained"
                if scr > 0
                else "accelerated"
                if scr < 0
                else "neutral"
            )
            results.append(
                {
                    "chrom": chrom,
                    "position": pos,
                    "phyloP_score": scr,
                    "status": status,
                }
            )

    return results


def distribute_blocks(blocks, n_processes):
    """
    Distribute blocks evenly across worker processes.

    Returns a list of block lists for each process.
    """
    assignments = [[] for _ in range(n_processes)]
    for i, block in enumerate(blocks):
        assignments[i % n_processes].append(block)
    return assignments


def apply_statistical_analysis(df, alpha=0.05, fdr_combined=False):
    """Apply FDR correction and significance labeling to phyloP scores."""
    logger.info(f"Applying FDR correction: alpha={alpha}, combined={fdr_combined}")

    scores = df["phyloP_score"].values
    p_values = np.power(10, -np.abs(scores))
    q_values = np.full_like(p_values, np.nan, dtype=float)
    significant = np.zeros_like(p_values, dtype=bool)

    if fdr_combined:
        _, q_all, _, _ = multipletests(p_values, alpha=alpha, method="fdr_bh")
        q_values[:] = q_all
        significant[:] = q_all <= alpha
    else:
        mask_pos = scores > 0
        mask_neg = scores < 0
        if mask_pos.any():
            _, q_pos, _, _ = multipletests(
                p_values[mask_pos], alpha=alpha, method="fdr_bh"
            )
            q_values[mask_pos] = q_pos
            significant[mask_pos] = q_pos <= alpha
        if mask_neg.any():
            _, q_neg, _, _ = multipletests(
                p_values[mask_neg], alpha=alpha, method="fdr_bh"
            )
            q_values[mask_neg] = q_neg
            significant[mask_neg] = q_neg <= alpha

    df["p_value"] = p_values
    df["q_value"] = q_values
    df["significant"] = np.where(significant, "yes", "no")
    return df


def generate_plots(df, output_file):
    """Generate diagnostic histograms and Q-Q plot for phyloP analysis."""
    plot_dir = os.path.dirname(output_file) or "."
    os.makedirs(plot_dir, exist_ok=True)
    prefix = os.path.splitext(os.path.basename(output_file))[0]

    scores = df["phyloP_score"].astype(float)
    pvals = df["p_value"].dropna().astype(float)

    # Histogram of scores
    plt.figure(figsize=(10, 6))
    scores.hist(bins=50, alpha=0.7)
    plt.axvline(0, color="r", linestyle="--")
    plt.title("Distribution of phyloP Scores")
    plt.xlabel("Score")
    plt.ylabel("Frequency")
    plt.savefig(f"{plot_dir}/{prefix}_score_histogram.png")
    plt.close()

    # Histogram of p-values
    plt.figure(figsize=(10, 6))
    pvals.hist(bins=20, alpha=0.7)
    plt.title("Distribution of p-values")
    plt.xlabel("p-value")
    plt.ylabel("Frequency")
    plt.savefig(f"{plot_dir}/{prefix}_pvalue_histogram.png")
    plt.close()

    # Q-Q plot
    sorted_p = np.sort(pvals)
    n = len(sorted_p)
    expected = np.linspace(0, 1, n + 2)[1:-1]
    plt.figure(figsize=(8, 8))
    plt.plot(expected, sorted_p, "o", markersize=2)
    plt.plot([0, 1], [0, 1], "r--")
    plt.title("Q-Q Plot of p-values")
    plt.xlabel("Expected p-value")
    plt.ylabel("Observed p-value")
    plt.savefig(f"{plot_dir}/{prefix}_pvalue_qq.png")
    plt.close()

    # Overlay histogram by significance
    plt.figure(figsize=(12, 6))
    sig = df[df["significant"] == "yes"]["phyloP_score"]
    not_sig = df[df["significant"] == "no"]["phyloP_score"]
    plt.hist([sig, not_sig], bins=50, alpha=0.7, label=["Significant", "Not significant"])
    plt.axvline(0, color="r", linestyle="--")
    plt.title("Score Distribution by Significance")
    plt.xlabel("phyloP Score")
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(f"{plot_dir}/{prefix}_score_by_significance.png")
    plt.close()

    logger.info(f"Diagnostic plots saved in {plot_dir}")


def print_summary(df):
    """Log and print summary statistics of the analysis."""
    total = len(df)
    sig = df[df["significant"] == "yes"]
    constr = sig[sig["status"] == "constrained"]
    accel = sig[sig["status"] == "accelerated"]

    logger.info("Analysis summary:")
    logger.info(f"Total positions: {total}")
    logger.info(f"Significant: {len(sig)} ({len(sig)/total*100:.2f}%)")
    logger.info(f"  - Constrained: {len(constr)} ({len(constr)/total*100:.2f}%)")
    logger.info(f"  - Accelerated: {len(accel)} ({len(accel)/total*100:.2f}%)")

    print("\nAnalysis summary:")
    print(f"Total positions: {total}")
    print(f"Significant: {len(sig)} ({len(sig)/total*100:.2f}%)")
    print(f"  - Constrained: {len(constr)} ({len(constr)/total*100:.2f}%)")
    print(f"  - Accelerated: {len(accel)} ({len(accel)/total*100:.2f}%)")


def parallel_process_phylop_file(
    input_file,
    output_file,
    alpha=0.05,
    plot=False,
    n_processes=None,
    fdr_combined=False,
):
    """
    Process phyloP scores file in parallel and optionally generate plots.

    Args:
        input_file: Path to phyloP WIG or TSV input.
        output_file: Path for tab-delimited output.
        alpha: Significance threshold for FDR correction.
        plot: If True, generate diagnostic plots.
        n_processes: Number of worker processes to use.
        fdr_combined: If True, apply FDR correction on all scores together.
    """
    start = time.time()
    n_processes = n_processes or min(mp.cpu_count(), 8)

    logger.info(
        f"Starting parallel processing: {input_file} -> {output_file} "
        f"using {n_processes} processes, FDR combined={fdr_combined}"
    )

    blocks, total_lines, data_lines = analyze_file_structure(input_file)
    assignments = distribute_blocks(blocks, n_processes)

    all_results = []
    if n_processes > 1:
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            futures = [
                executor.submit(process_blocks, input_file, chunk, idx)
                for idx, chunk in enumerate(assignments)
            ]
            for idx, fut in enumerate(futures):
                try:
                    all_results.extend(fut.result())
                except Exception as e:
                    logger.error(f"Error in chunk {idx}: {e}")
    else:
        for idx, chunk in enumerate(assignments):
            all_results.extend(process_blocks(input_file, chunk, idx))

    logger.info(
        f"Data extracted: {len(all_results)} positions (expected ~{data_lines})"
    )

    df = pd.DataFrame(all_results)
    if df.empty:
        logger.error("No valid data parsed; exiting.")
        return df

    df = apply_statistical_analysis(df, alpha=alpha, fdr_combined=fdr_combined)

    if plot:
        try:
            generate_plots(df, output_file)
        except Exception as e:
            logger.error(f"Plot generation failed: {e}")

    df.to_csv(output_file, sep="\t", index=False, float_format="%.6e")
    logger.info(f"Results saved to {output_file}")

    print_summary(df)

    elapsed = time.time() - start
    memory = psutil.Process(os.getpid()).memory_info().rss / 1024**2
    logger.info(f"Completed in {elapsed:.2f}s, peak memory {memory:.2f} MB")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Identify significant regions from phyloP scores."
    )
    parser.add_argument("input_file", help="phyloP input file path")
    parser.add_argument("output_file", help="output TSV file path")
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR significance level (default 0.05)",
    )
    parser.add_argument(
        "--plot", action="store_true", help="generate diagnostic plots"
    )
    parser.add_argument(
        "--processes",
        type=int,
        help="number of parallel processes (default=CPU count)",
    )
    parser.add_argument(
        "--fdr-combined",
        action="store_true",
        help="run FDR correction on entire dataset",
    )
    parser.add_argument("--debug", action="store_true", help="enable debug logging")

    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)

    parallel_process_phylop_file(
        args.input_file,
        args.output_file,
        alpha=args.alpha,
        plot=args.plot,
        n_processes=args.processes,
        fdr_combined=args.fdr_combined,
    )


if __name__ == "__main__":
    main()
