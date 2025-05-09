# Evolutionary constraint across the grass family (Poaceae)

## Background

The study of **conserved genomic elements** is a key direction in evolutionary biology. Genomic regions that have remained unchanged across millions of years often indicate essential biological functions. **Comparative genomics**, which involves comparing genomes across multiple species, provides a powerful framework to identify such regions under **purifying selection**, which eliminates deleterious mutations.

For example, the **Zoonomia Project**, a large-scale comparison of mammalian genomes, has demonstrated that this approach enables precise identification of conserved elements, including **non-coding regions**, that are likely to be functionally important. Identifying these constrained elements helps link genome structure to phenotype and function.

In the context of **grasses (Poaceae)** — and especially in species like **wheat (*Triticum aestivum*)** — this analysis is particularly relevant. The wheat genome is **large and TE-rich**, containing a vast number of **transposable elements (TEs)**. Despite this, there are **gene-rich regions** containing important loci, many of which have been under selection during domestication and breeding.

Studying conservation within these regions — and especially within **selective sweeps** — helps to:

- Understand **fundamental biological processes**,
- Reveal mechanisms of **adaptation** to different environmental and agricultural conditions,
- Identify genes responsible for key **agronomic traits**, such as:
  - **Brittle rachis** (shattering resistance),
  - **Free threshing** (grain harvestability).

Conservation of specific sequences in these functionally relevant regions, in contrast to rapidly evolving or highly repetitive regions, highlights their **critical role in plant biology** and crop improvement.


## Purpose

The purpose of this project is to explore evolutionary constraints and innovations across the grass family (Poaceae) using multispecies whole-genome alignments and phyloP conservation scores.

Students are provided with a reference-based multiple sequence alignment of 80 Poaceae genomes (reference: *Triticum aestivum*, bread wheat). The analysis is divided into the following goals:

- Identify regions of **conserved** and **accelerated** evolution across the genome using phyloP scores.
- Annotate conserved segments, most of which are expected to lie in **non-coding regions**.
- Investigate **gene-specific evolutionary rates** by comparing root-to-tip distances between species, and detect **clade-specific acceleration**, particularly between C₃ and C₄ lineages.
- Explore associations between evolutionary signals and **functional phenotypic traits**.

This is an exploratory project encouraging **independent hypothesis formulation** and **flexible workflows**. A high-performance server is available for conducting computationally intensive steps.

---

## Table of Contents

- [Description](##description)
- [Repository Structure](##repository-structure)
- [Results](##results)
- [Running the Analysis](##running-the-analysis)
- [Bibliography](##Bibliography)

---

## Description

This project includes two complementary analyses:

1. **Constraint analysis based on phyloP scores:**
   - Identifies 100kb genomic windows enriched for significantly conserved positions (phyloP q < 0.05).
   - Integrates annotations (CDS, exons, genes, TEs, ISBPs, GC content).
   - Produces Manhattan plots highlighting conserved and potentially regulatory regions.

2. **Root-to-tip evolutionary rate analysis (C3 vs C4):**
   - Constructs phylogenetic trees for each gene across species.
   - Compares root-to-tip distances between C3 and C4 species using Wilcoxon test.
   - Applies FDR correction to identify genes with significantly different rates.

---

## Repository Structure

### Environment Setup

Python

R packages


### Results

#### Constraint Analysis

**Output table:** `results/merged_bindata.100000.tsv`

| Column               | Description                                                             |
|----------------------|-------------------------------------------------------------------------|
| `chr`, `start`, `end`| 100kb genomic window coordinates                                        |
| `totalBases`         | Positions with any phyloP score                                         |
| `conservedBases`     | Bases in conserved segments (from q < 0.05 + score > 2.5)               |
| `cds.distinct.sum`   | Total CDS coverage length                                               |
| `exon.distinct.sum`  | Total exon coverage length                                              |
| `gene.count`         | Number of overlapping genes                                             |
| `intergenic.length`  | Bases outside gene annotations                                          |
| `GC`                 | GC content of the window                                                |
| `TE.count`           | Overlaps with transposable elements                                     |
| `ISBP.count`         | Overlaps with ISBP-class TE elements                                    |
| `totalConservedBases`| Significantly conserved positions (q < 0.05 & phyloP > 2.5)             |
| `meanConsScore`      | Mean phyloP for significant bases                                       |
| `gene`               | Comma-separated list of gene IDs                                        |

Manhattan plots show outlier windows (q < 0.05).  
Some constrained windows lack annotated genes (№47, `gene = NA`), indicating possible conserved non-coding elements.

![img](images/100kb_constraint_manhattan_q.png)

![img](images/100kb_constraint_manhattan_medianScore.png)

##### Plot Interpretation

| Plot                   | Description                                                                                      | Key Observations                                                                                           |
|------------------------|--------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------|
| `-log10(q-value)`      | Statistical significance of deviation from model prediction (dashed line = q-value threshold 0.05) | Most constrained genes (q-value < 0.05) are located on chromosomes Chr1A–4A. Other significant windows (ID 47) may represent conserved non-coding regions. |
| `medianConsScore`      | Robust median of phyloP values per 100kb window (measures typical conservation level)            | Constraint is distributed in wide blocks. High phyloP score ≠ always low q-value. Top 10 constrained windows per chromosome are labeled. |

---

#### Evolutionary Analysis of Genes Associated with C₃/C₄ Photosynthesis


C₃ photosynthesis is the most common carbon fixation pathway, especially in temperate climates. It relies on the enzyme **RuBisCO**, which is prone to **photorespiration** under high temperatures or low CO₂ concentrations.

C₄ photosynthesis is an adaptation to **hot and dry environments**. It begins with the fixation of CO₂ into 4-carbon compounds in mesophyll cells. These compounds are transported to **bundle-sheath cells**, where CO₂ is released into a high-concentration zone for RuBisCO. This **reduces photorespiration** and increases efficiency.

##### Summary of Evolutionary Analysis

We analyzed evolutionary divergence between C₃ and C₄ species by calculating **root-to-tip distances (RTT)** in gene trees constructed from aligned sequences.  
Statistical comparison was performed using **Wilcoxon tests**, with **FDR correction** applied to obtain q-values.

**Results:**

- **220 genes** showed significant RTT divergence between C₃ and C₄ lineages (**q < 0.05**)
- **173 genes** had **longer RTT in C₃ species**  
  → Indicative of **faster evolution** in C₃ lineages
- **47 genes** had **longer RTT in C₄ species**  
  → Potentially **conserved** in C₄ lineages
- **No genes** showed opposite-direction significance

#### Functional Enrichment of C₃-Shifted Genes

C₃-associated genes (i.e., evolving faster in C₃) were significantly enriched for the following functions:

- **Oxidoreductase activity**  
  Involved in redox metabolism, key for photosynthetic electron transport
- **DNA binding and transcription regulation**  
  Includes transcription factors and DNA-interacting proteins
- **Iron and heme binding**  
  Important for cytochromes and electron transport components
- **Oxidation–reduction processes**  
  General redox-related cellular pathways
- **Protein interactions**  
  Involved in signal transduction and complex formation
- **Hydrolase activity**  
  Enzymes that hydrolyze chemical bonds

These categories reflect core **photosynthetic and regulatory functions**. Their **evolutionary constraint in C₄ species** may result from the compartmentalized and efficient operation of photosynthesis under **arid and high-light** conditions.


**GO enrichment** for C3-shifted genes highlights biological categories such as:

- Oxidoreductase activity  
- DNA binding  
- Transcription factor activity  
- Heme/iron ion binding

![img](images/boxplots_top_genes.png)

![img](images/Category_Heatmap_C3_C4.png)

#### Biological Interpretation

- **Genes with longer RTT in C₄** species may represent **adaptations** that evolved under **increased constraint** to maintain optimized function.
- **Genes with longer RTT in C₃** species likely underwent **relaxed selection** or **adaptive divergence** in C₃ lineages.

##### 1. Process phyloP data

```bash
python scripts/process_phyloP_paral.py data/roast.maf
```
This script processes phyloP scores from a MAF file and applies FDR correction to identify significantly conserved positions (q < 0.05)

##### 2. Define conserved regions

```bash
python scripts/find_conserved_multi.py
```
Merges adjacent significant positions into conserved segments. Output: {CHR}_conserved_segments.bed.

##### 3. Generate 100kb bins with annotations

```bash
bash scripts/generate_bindata.sh
```
This script:

Creates 100kb windows using bedtools makewindows

Computes the following features per bin:

Total number of bases with phyloP scores

Number of conserved positions

Length of overlapping CDS and exons

Gene count

Intergenic length

GC content (via bedtools nuc)

Number of overlaps with TE and ISBP elements

Mean phyloP score for significant bases

List of overlapping gene names

Output: merged_bindata.100000.tsv

##### 4. Identify constrained windows

```bash
Rscript scripts/100kb_bin_constraint.R
```
Fits a linear regression model to predict the expected number of conserved bases in each window, then computes standardized residuals, p-values, and q-values. Windows with q < 0.05 are considered significantly constrained. A Manhattan plot is generated to visualize results.

#### C3/C4 Root-to-Tip Analysis

##### 1. Prepare per-gene alignments

Use gene coordinates and the MAF file to extract one alignment per gene. Convert each MAF block into FASTA format using

```bash
mafToFastaStitcher input.gene.maf > gene.fasta
```
Ensure that species identifiers match between alignment files and the tree.

##### 2. Run evolutionary rate comparison

```bash
Rscript scripts/run_C3_C4.R | tee run_C3_C4.log
```
This script:

Builds gene trees using the JC69 model  
Computes root-to-tip (RTT) distances per species  
Compares RTT between C3 and C4 species using the Wilcoxon test  
Applies FDR correction to identify significantly shifted genes  

**Outputs:**

- `genes_C3.txt`: genes evolving faster in C3 lineages  
- `genes_C4.txt`: genes evolving faster in C4 lineages  
- `root_to_tip_C3vsC4.csv`: full table with RTT values, p-values, q-values  


## Bibliography

1. **Zoonomia Consortium**. (2020). A comparative genomics multitool for scientific discovery and conservation. *Nature*, 587, 240–245. https://doi.org/10.1038/s41586-020-2876-6

2. **Avraham A. Levy**, **Moshe Feldman**. (2022). Evolution and origin of bread wheat. *The Plant Cell*, 34, 2549–2567. https://doi.org/10.1093/plcell/koac130

3. **Hong Cheng**, **Jing Liu**, **Jia Wen**, **Xiaojun Nie**, et al. (2019). Frequent intra- and inter-species introgression shapes the landscape of genetic variation in bread wheat. *BMC Genomics*. https://doi.org/10.1186/s13059-019-1744-x

4. **Matthew J. Christmas**, **E. C. Hare**, **T. L. Capellini**, et al. (2023). Evolutionary constraint and innovation across hundreds of placental mammals. *Science*, 380, eabn3943. https://doi.org/10.1126/science.abn3943

