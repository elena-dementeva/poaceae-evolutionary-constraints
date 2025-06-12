import argparse
import os
import re
import subprocess
import tarfile
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import pandas as pd
from intervaltree import Interval, IntervalTree


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-g", "--gff", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("-t", "--threshold", type=float, default=1.5)
    parser.add_argument("--min-q", type=float, default=0.05)
    parser.add_argument("--max-gap", type=int, default=1000)
    parser.add_argument("--max-bad", type=int, default=0)
    parser.add_argument("--out", default=None)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=500_000)
    parser.add_argument("--overlap", type=int, default=10_000)
    parser.add_argument("--chunk-threads", type=int, default=cpu_count())
    parser.add_argument("--summary-threads", type=int, default=4)
    return parser.parse_args()


def load_phyloP(filename, chromosome, min_q_value):
    print(f"Loading & filtering {filename} for {chromosome} in chunks…")
    if filename.endswith((".tar.gz", ".tgz")):
        archive = tarfile.open(filename, "r:gz")
        member = next(
            (m for m in archive.getmembers() if m.name.endswith(".tsv")),
            None,
        )
        if member is None:
            raise ValueError(f"No .tsv inside {filename}")
        header_file = archive.extractfile(member)
    else:
        header_file = filename

    header = pd.read_csv(header_file, sep="\t", nrows=0).columns
    header_map = {c.lower(): c for c in header}

    chrom_col = next(
        (header_map[c] for c in ("chrom", "chr", "chromosome") if c in header_map),
        None,
    )
    pos_col = next(
        (header_map[c] for c in ("position", "pos") if c in header_map),
        None,
    )
    score_col = next(
        (c for c in header if re.search(r"phyloP", c, re.IGNORECASE)),
        None,
    )
    q_col = next(
        (c for c in header if re.search(r"\bq(?:_?value)?\b", c, re.IGNORECASE)),
        None,
    )

    if not pos_col or not score_col or not q_col:
        raise ValueError(
            f"Required columns not found: pos={pos_col}, score={score_col}, q={q_col}"
        )

    use_columns = [col for col in (chrom_col, pos_col, score_col, q_col) if col]
    dtypes = {
        pos_col: int,
        score_col: float,
        q_col: float,
    }

    if filename.endswith((".tar.gz", ".tgz")):
        source = archive.extractfile(member)
    else:
        source = filename

    reader = pd.read_csv(
        source,
        sep="\t",
        usecols=use_columns,
        dtype=dtypes,
        chunksize=500_000,
        low_memory=True,
    )

    filtered = []
    total_kept = 0
    for chunk in reader:
        if chrom_col and chromosome in chunk[chrom_col].unique():
            chunk = chunk[chunk[chrom_col] == chromosome]
        else:
            chunk["chrom"] = chromosome

        sel = chunk[(chunk[q_col] <= min_q_value) & (chunk[score_col] > 0)]
        total_kept += len(sel)
        if not sel.empty:
            sel = sel.rename(
                columns={
                    chrom_col or pos_col: "chrom",
                    pos_col: "position",
                    score_col: "phyloP_-log_pvalue",
                    q_col: "q_value",
                }
            )
            filtered.append(sel[["chrom", "position", "phyloP_-log_pvalue", "q_value"]])
        print(f"  chunk done, kept {len(sel)} rows")

    if not filtered:
        raise ValueError(
            "No positions kept after filtering; check column names or --chrom/--min-q"
        )

    df = pd.concat(filtered, ignore_index=True)
    print(f"Finished: {total_kept} positions after filter")
    return df


def make_chunks(df, chunk_size, overlap):
    for start_idx in range(0, len(df), chunk_size):
        lo = max(0, start_idx - overlap)
        hi = min(len(df), start_idx + chunk_size + overlap)
        yield df.iloc[lo:hi].reset_index(drop=True)


def find_continuous_segments(df_chunk, threshold, min_q, max_gap, max_bad):
    df_chunk = df_chunk.sort_values("position").reset_index(drop=True)
    positions = df_chunk["position"].to_numpy()
    scores = df_chunk["phyloP_-log_pvalue"].to_numpy()
    q_values = df_chunk["q_value"].to_numpy()

    segments = []
    idx = 0
    while idx < len(df_chunk):
        if q_values[idx] <= min_q and scores[idx] > 0:
            start_idx = idx
            end_idx = idx
            total_score = scores[idx]
            valid_count = 1
            bad_counter = 0
            while (
                end_idx + 1 < len(df_chunk)
                and positions[end_idx + 1] - positions[end_idx] <= max_gap
            ):
                nq = q_values[end_idx + 1]
                ns = scores[end_idx + 1]
                if nq > min_q or ns <= 0:
                    bad_counter += 1
                    if bad_counter > max_bad:
                        break
                else:
                    total_score += ns
                    valid_count += 1
                    bad_counter = 0
                end_idx += 1
            segments.append(
                (
                    df_chunk.iloc[0]["chrom"],
                    positions[start_idx] - 1,
                    positions[end_idx],
                    float(total_score / valid_count),
                )
            )
            idx = end_idx + 1
        else:
            idx += 1
    print(f"   chunk: found {len(segments)} segments")
    return pd.DataFrame(segments, columns=["chrom", "start", "end", "mean_score"])


def load_and_sort_beds(gff_file, prefix, chromosome):
    features = {}
    print("Parsing GFF and saving feature BED files")
    for feature in ("CDS", "exon", "gene"):
        records = []
        with open(gff_file) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip().split("\t")
                if parts[0] == chromosome and parts[2] == feature:
                    records.append((parts[0], int(parts[3]) - 1, int(parts[4])))
        df_feat = pd.DataFrame(records, columns=["chrom", "start", "end"])
        bed_file = f"{prefix}_{feature}.bed"
        df_feat.to_csv(bed_file, sep="\t", header=False, index=False)
        subprocess.run(f"sort -k1,1 -k2,2n {bed_file} -o {bed_file}", shell=True)
        print(f"   {feature}: {len(df_feat)} intervals → {bed_file} (sorted)")
        features[feature] = df_feat
    return features


def summarize_feature(feature, segs_df, feats_df, prefix, no_plots):
    print(f"\nFeature: {feature}")
    print("Summary:")
    print(f"   All segments: {len(segs_df)}")
    if feature != "noncoding":
        tree = IntervalTree(
            Interval(r.start, r.end + 1) for _, r in feats_df.iterrows()
        )
        hits = sum(
            bool(tree.overlap(r.start, r.end + 1)) for _, r in segs_df.iterrows()
        )
        pct = hits / len(segs_df) * 100 if len(segs_df) else 0
        print(f"   Intersected with {feature}: {hits} ({pct:.1f}%)")
    lengths = segs_df.end - segs_df.start + 1
    q25, q50, q75, q95 = lengths.quantile([0.25, 0.5, 0.75, 0.95])
    print(
        f"   Length (bp): 25%≤{int(q25)}, med={int(q50)}, 75%≤{int(q75)}, "
        f"95%≤{int(q95)}"
    )
    if not no_plots:
        mids = ((segs_df.start + segs_df.end) // 2) / 1_000_000
        plt.figure(figsize=(10, 4))
        plt.scatter(
            mids,
            segs_df.mean_score,
            s=12,
            alpha=0.6,
            edgecolor="none",
            color="steelblue",
        )
        plt.title(f"{feature} Manhattan – {segs_df.iloc[0]['chrom']}")
        plt.xlabel("Genomic position (Mb)")
        plt.ylabel("Mean phyloP")
        plt.ylim(0, 5)
        for x in range(int(mids.max()) + 1):
            plt.axvline(x=x, linestyle="--", alpha=0.2)
        plt.tight_layout()
        plt.savefig(f"{prefix}_{feature}_manhattan.png", dpi=300)
        plt.close()


def main():
    args = parse_args()
    prefix = args.out or os.path.splitext(os.path.basename(args.input))[0]

    df_phy = load_phyloP(args.input, args.chrom, args.min_q)

    chunks = list(make_chunks(df_phy, args.chunk_size, args.overlap))
    print(f"Split input into {len(chunks)} chunks")
    print(f"Processing chunks with {args.chunk_threads} threads")

    with Pool(args.chunk_threads) as pool:
        dfs = pool.starmap(
            find_continuous_segments,
            [
                (
                    chunk,
                    args.threshold,
                    args.min_q,
                    args.max_gap,
                    args.max_bad,
                )
                for chunk in chunks
            ],
        )

    print("Segments processed in parallel")
    segs_df = pd.concat(dfs, ignore_index=True).drop_duplicates(
        ["chrom", "start", "end"]
    )
    print(f"Unique segments: {len(segs_df)}")

    print("Saving segment tables")
    segs_tsv = f"{prefix}_segments.tsv"
    segs_df.to_csv(segs_tsv, sep="\t", index=False)
    seg_bed = f"{prefix}_segments.bed"
    segs_df[["chrom", "start", "end"]].to_csv(
        seg_bed, sep="\t", header=False, index=False
    )
    subprocess.run(f"sort -k1,1 -k2,2n {seg_bed} -o {seg_bed}", shell=True)

    feats = load_and_sort_beds(args.gff, prefix, args.chrom)

    print("Calculating noncoding regions and exporting TSV")
    non_bed = f"{prefix}_noncoding.bed"
    beds = " ".join(f"-b {prefix}_{feat}.bed" for feat in feats)
    subprocess.run(f"bedtools subtract -A -a {seg_bed} {beds} > {non_bed}", shell=True)

    non_df = pd.read_csv(
        non_bed, sep="\t", header=None, names=["chrom", "start", "end"]
    )
    noncoding_df = pd.merge(non_df, segs_df, on=["chrom", "start", "end"])
    non_tsv = f"{prefix}_noncoding_segments.tsv"
    noncoding_df.to_csv(non_tsv, sep="\t", index=False)
    print(f"   Noncoding TSV → {non_tsv}")

    for feat in ["CDS", "exon", "gene"]:
        summarize_feature(feat, segs_df, feats[feat], prefix, args.no_plots)
    summarize_feature("noncoding", noncoding_df, None, prefix, args.no_plots)

    print("\nDone.")


if __name__ == "__main__":
    main()
