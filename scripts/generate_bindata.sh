#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -g GENOME_FA -f GFF -t TE_GFF -i ISBP_BED -c ANNOT_CSV --tsv1 CHR1A_TSV ... --tsv7 CHR7A_TSV"
  exit 1
}

TSV_FILES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g) GENOME_FA=$2; shift 2 ;;
    -f) GFF=$2; shift 2 ;;
    -t) TE_GFF=$2; shift 2 ;;
    -i) ISBP_BED=$2; shift 2 ;;
    -c) ANNOT_CSV=$2; shift 2 ;;
    --tsv1) TSV_FILES[0]=$2; shift 2 ;;
    --tsv2) TSV_FILES[1]=$2; shift 2 ;;
    --tsv3) TSV_FILES[2]=$2; shift 2 ;;
    --tsv4) TSV_FILES[3]=$2; shift 2 ;;
    --tsv5) TSV_FILES[4]=$2; shift 2 ;;
    --tsv6) TSV_FILES[5]=$2; shift 2 ;;
    --tsv7) TSV_FILES[6]=$2; shift 2 ;;
    *) usage ;;
  esac
done

# Check required inputs
: "${GENOME_FA:?Missing -g}"
: "${GFF:?Missing -f}"
: "${TE_GFF:?Missing -t}"
: "${ISBP_BED:?Missing -i}"
: "${ANNOT_CSV:?Missing -c}"
for i in {0..6}; do
  if [[ -z "${TSV_FILES[$i]:-}" ]]; then
    echo "Missing --tsv$((i+1))"
    exit 1
  fi
done

ANNOT_TSV=annotations.tsv
GENOME_SIZE=TriticumA.genome
SCRIPT_DIR=$(pwd)
CHR_LIST=(Chr1A Chr2A Chr3A Chr4A Chr5A Chr6A Chr7A)

# Prepare annotation table
if [[ ! -f $ANNOT_TSV ]]; then
  echo "Preparing annotation table (gene_id, name, description, GO)..."
  awk -F',' 'NR>1 {
    id = $1
    if ($2 ~ /Pfam|InterPro|Human readable description/) {
      if (name[id] == "" && $5 != "") name[id] = $5
    }
    if ($2 == "Human readable description" && desc[id] == "" && $4 != "") {
      desc[id] = $4
    }
    if ($2 == "Gene Ontology" && $5 != "") {
      go[id] = go[id] "," $5
    }
  }
  END {
    for (id in name) {
      print id "\t" name[id] "\t" desc[id] "\t" substr(go[id],2)
    }
  }' "$ANNOT_CSV" | sort > "$ANNOT_TSV"
else
  echo "Annotation file already exists: using $ANNOT_TSV"
fi

# Check genome index
echo "Checking if genome index exists..."
if [[ ! -f "${GENOME_FA}.fai" ]]; then
  echo "Index not found. Indexing genome..."
  samtools faidx "$GENOME_FA"
fi

cut -f1,2 "${GENOME_FA}.fai" > "$GENOME_SIZE"

# Loop over chromosomes
for i in "${!CHR_LIST[@]}"; do
CHR="${CHR_LIST[$i]}"
TSV="${TSV_FILES[$i]}"
(
  echo "======== $CHR: START ========"

  if [[ ! -f "$TSV" ]]; then
    echo "WARNING: $TSV not found, skipping $CHR"
    exit 0
  fi

  OUTDIR="${CHR}_bindata"
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"

  echo "Generating bins and copying conserved segments for $CHR"
  bedtools makewindows -g "$SCRIPT_DIR/$GENOME_SIZE" -w 100000 > bins100k.bed
  cp "$SCRIPT_DIR/${CHR}_conserved/${CHR}_conserved_segments.bed" cons.bed

  echo "Filtering all phyloP and conserved phyloP positions for $CHR"
  awk 'NR>1 {print $1"\t"$2-1"\t"$2}' "$TSV" > all_phyloP.bed
  awk 'NR>1 && $6 < 0.05 && $3 > 2.5 {print $1"\t"($2-1)"\t"$2"\t"$3}' "$TSV" > cons_signif_phyloP.bed

  echo "Extracting annotation features for $CHR"
  awk '$3 == "CDS" {print $1"\t"$4-1"\t"$5}' "$GFF" > CDS.bed
  awk '$3 == "exon" {print $1"\t"$4-1"\t"$5}' "$GFF" > exon.bed
  awk -F'\t' '$3 == "gene" {
    match($9, /ID=([^;]+)/, id);
    name = "NA";
    if ($9 ~ /Name=/) {
      match($9, /Name=([^;]+)/, n);
      name = n[1];
    }
    print $1 "\t" $4-1 "\t" $5 "\t" name;
  }' "$GFF" > gene_with_names.bed
  awk '$3 == "transposable_element" {print $1"\t"$4-1"\t"$5}' "$TE_GFF" > TE.bed
  cp "$ISBP_BED" ISBP.bed

  echo "Sorting files for $CHR"
  for file in all_phyloP cons_signif_phyloP bins100k cons CDS exon gene_with_names TE ISBP; do
    sort -k1,1 -k2,2n "$file.bed" > "$file.sorted.bed"
  done

  echo "Calculating GC content for $CHR"
  bedtools nuc -fi "$GENOME_FA" -bed bins100k.bed > gc_stats.txt

  echo "Computing coverage and intersections for $CHR"
  bedtools coverage -a bins100k.sorted.bed -b all_phyloP.sorted.bed -sorted -counts > all.cov
  bedtools coverage -a bins100k.sorted.bed -b cons.sorted.bed -sorted -counts > cons.cov
  bedtools coverage -a bins100k.sorted.bed -b CDS.sorted.bed -sorted -counts > cds.cov
  bedtools coverage -a bins100k.sorted.bed -b exon.sorted.bed -sorted -counts > exon.cov
  bedtools intersect -a bins100k.sorted.bed -b gene_with_names.sorted.bed -sorted -c > gene.count.cov
  bedtools coverage -a bins100k.sorted.bed -b TE.sorted.bed -sorted -counts > TE.cov
  bedtools coverage -a bins100k.sorted.bed -b ISBP.sorted.bed -sorted -counts > ISBP.cov
  bedtools subtract -a bins100k.sorted.bed -b gene_with_names.sorted.bed -sorted > intergenic.tmp.bed
  bedtools coverage -a bins100k.sorted.bed -b intergenic.tmp.bed -sorted -counts > intergenic.cov

  echo "Calculating totalConservedBases and meanConsScore for $CHR"
  bedtools coverage -a bins100k.sorted.bed -b cons_signif_phyloP.sorted.bed -sorted -counts > cons_signif.count.cov
  bedtools map -a bins100k.sorted.bed -b cons_signif_phyloP.sorted.bed -c 4 -o mean > cons_signif.mean.bed

  echo "Generating gene field with annotations for $CHR"
  bedtools intersect -a bins100k.sorted.bed -b gene_with_names.sorted.bed -wa -wb |
    awk '{print $1"\t"$2"\t"$3"\t"$NF}' | sort -k1,1 -k2,2n > gene.names.bed

  awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3; genes[key]=genes[key]","$4}
       END{for (k in genes) {print k"\t"substr(genes[k], 2)}}' gene.names.bed |
       sort -k1,1 -k2,2n > gene_field.tsv

  echo "Annotating gene field for $CHR"
  awk 'NR==FNR {
    ann[$1] = $0
    next
  }
  {
    split($0, f, "\t")
    id = f[4]
    if (id in ann) {
      split(ann[id], a, "\t")
      name = a[2]; desc = a[3]; go = a[4]
    } else {
      name = id; desc = "NA"; go = "NA"
    }
    print $0 "\t" name "\t" desc "\t" go
  }' "$SCRIPT_DIR/$ANNOT_TSV" gene_field.tsv > gene_field.annot.tsv

  echo "Writing final bindata.100000.tsv for $CHR"
  paste \
    <(cut -f1-3 bins100k.sorted.bed) \
    <(cut -f4 all.cov) <(cut -f4 cons.cov) <(cut -f4 cds.cov) <(cut -f4 exon.cov) \
    <(cut -f4 gene.count.cov) <(cut -f4 intergenic.cov) <(tail -n +2 gc_stats.txt | cut -f5) \
    <(cut -f4 TE.cov) <(cut -f4 ISBP.cov) <(cut -f4 cons_signif.count.cov) <(cut -f4 cons_signif.mean.bed) \
    <(cut -f4 gene_field.annot.tsv) <(cut -f5 gene_field.annot.tsv) <(cut -f6 gene_field.annot.tsv) <(cut -f7 gene_field.annot.tsv) \
    > bindata.100000.tsv

  echo "======== $CHR: FINISHED ========"
) &
done

wait

echo "Merging all *_bindata/bindata.100000.tsv files..."
first=1
> merged_bindata.100000.tsv
for f in *_bindata/bindata.100000.tsv; do
  if [ $first -eq 1 ]; then
    head -n 1 "$f" > merged_bindata.100000.tsv
    first=0
  fi
  tail -n +2 "$f" >> merged_bindata.100000.tsv
done

echo "Done! Output: merged_bindata.100000.tsv"
