import argparse
import os
import re
import subprocess
import tarfile
import pandas as pd
import matplotlib.pyplot as plt
from intervaltree import IntervalTree, Interval
from multiprocessing import Pool, cpu_count


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-g", "--gff", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("-t", "--threshold", type=float, default=1.5)
    p.add_argument("--min-q", type=float, default=0.05)
    p.add_argument("--max-gap", type=int, default=1000)
    p.add_argument("--max-bad", type=int, default=0)
    p.add_argument("--out", default=None)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--chunk-size", type=int, default=500_000)
    p.add_argument("--overlap", type=int, default=10_000)
    p.add_argument("--chunk-threads", type=int, default=cpu_count())
    p.add_argument("--summary-threads", type=int, default=4)
    return p.parse_args()


def load_phyloP(fn, chrom, min_q):
    print(f"Loading & filtering {fn} for {chrom} in chunks…")
    if fn.endswith((".tar.gz", ".tgz")):
        t = tarfile.open(fn, "r:gz")
        member = next((m for m in t.getmembers() if m.name.endswith(".tsv")), None)
        if member is None:
            raise ValueError(f"No .tsv inside {fn}")
        header_f = t.extractfile(member)
    else:
        header_f = fn
    header = pd.read_csv(header_f, sep="\t", nrows=0).columns.tolist()

    chrom_col = next((c for c in header if c.lower() in ("chrom","chr","chromosome")), None)
    pos_col   = next((c for c in header if c.lower() in ("position","pos")), None)
    score_col = next((c for c in header if re.match(r"phyloP", c, re.IGNORECASE)), None)
    q_col     = next((c for c in header if "q" in c.lower()), None)

    if not pos_col or not score_col or not q_col:
        raise ValueError(f"Required columns not found: pos={pos_col}, score={score_col}, q={q_col}")

    usecols = [col for col in (chrom_col, pos_col, score_col, q_col) if col]
    dtypes  = {pos_col: int, score_col: float, q_col: float}

    if fn.endswith((".tar.gz", ".tgz")):
        source = t.extractfile(member)
    else:
        source = fn

    reader = pd.read_csv(
        source,
        sep="\t",
        usecols=usecols,
        dtype=dtypes,
        chunksize=500_000,
        low_memory=True
    )

    kept = []
    total = 0
    for chunk in reader:
        if chrom_col and chrom in chunk[chrom_col].unique():
            chunk = chunk[chunk[chrom_col] == chrom]
        else:
            chunk["chrom"] = chrom

        sel = chunk[(chunk[q_col] <= min_q) & (chunk[score_col] > 0)]
        total += len(sel)
        if not sel.empty:
            sel = sel.rename(columns={
                chrom_col or pos_col: "chrom",
                pos_col:   "position",
                score_col: "phyloP_-log_pvalue",
                q_col:     "q_value"
            })
            kept.append(sel[["chrom","position","phyloP_-log_pvalue","q_value"]])
        print(f"  chunk done, kept {len(sel)} rows")

    if not kept:
        raise ValueError("No positions kept after filtering — проверьте имена колонок или параметры `--chrom`/`--min-q`")

    df = pd.concat(kept, ignore_index=True)
    print(f"Finished: {total} positions after filter")
    return df


def make_chunks(df, size, overlap):
    for start in range(0, len(df), size):
        lo = max(0, start - overlap)
        hi = min(len(df), start + size + overlap)
        yield df.iloc[lo:hi].reset_index(drop=True)

def find_continuous_segments(df_chunk, thr, min_q, max_gap, max_bad):
    df_chunk = df_chunk.sort_values("position").reset_index(drop=True)
    pos = df_chunk["position"].to_numpy()
    score = df_chunk["phyloP_-log_pvalue"].to_numpy()
    qv = df_chunk["q_value"].to_numpy()
    regs = []
    i = 0
    while i < len(df_chunk):
        if qv[i] <= min_q and score[i] > 0:
            s, e = i, i
            tot, cnt = score[i], 1
            bad = 0
            while e + 1 < len(df_chunk) and pos[e + 1] - pos[e] <= max_gap:
                new_q, new_s = qv[e + 1], score[e + 1]
                if new_q > min_q or new_s <= 0:
                    bad += 1
                    if bad > max_bad:
                        break
                else:
                    tot += new_s
                    cnt += 1
                    bad = 0
                e += 1
            regs.append((df_chunk.iloc[0]["chrom"], pos[s] - 1, pos[e], float(tot / cnt)))
            i = e + 1
        else:
            i += 1
    print(f"   chunk: found {len(regs)} segments")
    return pd.DataFrame(regs, columns=["chrom", "start", "end", "mean_score"])

def load_and_sort_beds(gff, prefix, chrom):
    feats = {}
    print("6) Parsing GFF and saving feature BED files")
    for feat in ("CDS", "exon", "gene"):
        rec = []
        with open(gff) as f:
            for L in f:
                if L.startswith("#"): continue
                parts = L.rstrip().split("\t")
                if parts[0] == chrom and parts[2] == feat:
                    rec.append((parts[0], int(parts[3]) - 1, int(parts[4])))
        df = pd.DataFrame(rec, columns=["chrom","start","end"])
        bed = f"{prefix}_{feat}.bed"
        df.to_csv(bed, sep="\t", header=False, index=False)
        subprocess.run(f"sort -k1,1 -k2,2n {bed} -o {bed}", shell=True)
        print(f"   {feat}: {len(df)} intervals → {bed} (sorted)")
        feats[feat] = df
    return feats

def summarize_feature(feature, segs_df, feats_df, prefix, no_plots):
    print(f"\nfeature: {feature}")
    print("7) Sunnary:")
    print(f"   All segments: {len(segs_df)}")
    if feature != "noncoding":
        tree = IntervalTree(Interval(r.start, r.end+1) for _, r in feats_df.iterrows())
        hits = sum(bool(tree.overlap(r.start, r.end+1)) for _, r in segs_df.iterrows())
        pct = hits / len(segs_df) * 100 if len(segs_df) else 0
        print(f"   Interlined with {feature}: {hits} ({pct:.1f}%)")
    lengths = segs_df.end - segs_df.start + 1
    q25, q50, q75, q95 = lengths.quantile([0.25, 0.5, 0.75, 0.95])
    print(f"   Lengh (bp): 25%≤{int(q25)}, med={int(q50)}, 75%≤{int(q75)}, 95%≤{int(q95)}")
    if not no_plots:
        mids = ((segs_df.start + segs_df.end) // 2) / 1_000_000
        plt.figure(figsize=(10,4))
        plt.scatter(mids, segs_df.mean_score, s=12, alpha=0.6, edgecolor='none', color='steelblue')
        plt.title(f"{feature} Manhattan – {segs_df.iloc[0]['chrom']}")
        plt.xlabel("Genomic position (Mb)")
        plt.ylabel("Mean phyloP")
        plt.ylim(0, 5)
        for x in range(0, int(mids.max()) + 1):
            plt.axvline(x=x, color='gray', linestyle='--', alpha=0.2)
        plt.tight_layout()
        plt.savefig(f"{prefix}_{feature}_manhattan.png", dpi=300)
        plt.close()

def main():
    args = parse_args()
    prefix = args.out or os.path.splitext(os.path.basename(args.input))[0]

    df_phy = load_phyloP(args.input, args.chrom, args.min_q)

    chunks = list(make_chunks(df_phy, args.chunk_size, args.overlap))
    print(f"2) Split input into {len(chunks)} chunks")
    print(f"3) Processing chunks with {args.chunk_threads} threads")

    with Pool(args.chunk_threads) as pool:
        dfs = pool.starmap(find_continuous_segments, [(ch, args.threshold, args.min_q, args.max_gap, args.max_bad) for ch in chunks])

    print("4) Segments processed in parallel")
    segs_df = pd.concat(dfs, ignore_index=True).drop_duplicates(["chrom", "start", "end"])
    print(f"   Unique segments: {len(segs_df)}")

    print("5) Saving segment tables")
    segs_tsv = f"{prefix}_segments.tsv"
    segs_df.to_csv(segs_tsv, sep="\t", index=False)
    seg_bed = f"{prefix}_segments.bed"
    segs_df[["chrom","start","end"]].to_csv(seg_bed, sep="\t", header=False, index=False)
    subprocess.run(f"sort -k1,1 -k2,2n {seg_bed} -o {seg_bed}", shell=True)

    feats = load_and_sort_beds(args.gff, prefix, args.chrom)

    print("7) Calculating noncoding regions and exporting TSV")
    non_bed = f"{prefix}_noncoding.bed"
    beds = " ".join(f"-b {prefix}_{feat}.bed" for feat in feats)
    subprocess.run(f"bedtools subtract -A -a {seg_bed} {beds} > {non_bed}", shell=True)

    non_df = pd.read_csv(non_bed, sep="\t", header=None, names=["chrom","start","end"])
    noncoding_df = pd.merge(non_df, segs_df, on=["chrom","start","end"])
    non_tsv = f"{prefix}_noncoding_segments.tsv"
    noncoding_df.to_csv(non_tsv, sep="\t", index=False)
    print(f"   Noncoding TSV → {non_tsv}")

    for feat in ["CDS", "exon", "gene"]:
        summarize_feature(feat, segs_df, feats[feat], prefix, args.no_plots)
    summarize_feature("noncoding", noncoding_df, None, prefix, args.no_plots)
    print("\nDone.")

if __name__ == "__main__":
    main()
