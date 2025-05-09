library(ape)
library(qvalue)
library(ggplot2)
library(reshape2)
library(phytools)

clean_names <- function(x) sub("\\..*$", "", trimws(x))

cat("Loading phenotype data...\n")
pheno <- read.csv("species_c4_table.csv", colClasses = c("character", "numeric"))
pheno[[1]] <- clean_names(pheno[[1]])
pheno_vec  <- setNames(pheno[[2]], pheno[[1]])

cat("Loading FASTA files and phylogenetic tree...\n")
fasta_files <- list.files("gene_fastas", "\\.fa(sta)?$", full.names = TRUE)
tree <- read.tree("tree_poaceae.nwk")

pvals <- ntips <- numeric()
rtt_diff <- numeric()

cat("Analyzing genes...\n")
for (f in fasta_files) {
  gene <- basename(f)
  dna  <- read.dna(f, "fasta")
  rownames(dna) <- clean_names(rownames(dna))
  spp <- intersect(rownames(dna), names(pheno_vec))
  if (length(spp) < 4) next
  
  dmat <- dist.dna(dna[spp, ], "JC69", pairwise.deletion = TRUE)
  if (any(is.na(dmat))) next
  
  gtree <- njs(dmat)
  rtt <- node.depth.edgelength(gtree)[match(gtree$tip.label, gtree$tip.label)]
  names(rtt) <- gtree$tip.label
  
  c3 <- rtt[names(pheno_vec[pheno_vec == 0])]
  c4 <- rtt[names(pheno_vec[pheno_vec == 1])]
  c3 <- c3[!is.na(c3)]
  c4 <- c4[!is.na(c4)]
  if (length(c3) < 2 || length(c4) < 2) next
  
  pvals[gene] <- wilcox.test(c3, c4)$p.value
  ntips[gene] <- length(rtt)
  rtt_diff[gene] <- mean(c4, na.rm = TRUE) - mean(c3, na.rm = TRUE)
}

stopifnot(length(pvals) > 0)
qvals <- qvalue(pvals)$qvalues

cat("\nTop p-values:\n")
print(head(sort(pvals), 10))

cat("\nTop q-values after FDR correction:\n")
print(head(sort(qvals), 10))

cat("\nRunning phylANOVA on the full tree...\n")
common <- intersect(tree$tip.label, names(pheno_vec))
tree_pruned <- keep.tip(tree, common)
pheno_vec_filtered <- pheno_vec[common]

rtt_all <- node.depth.edgelength(tree_pruned)[match(tree_pruned$tip.label, tree_pruned$tip.label)]
names(rtt_all) <- tree_pruned$tip.label
groups <- factor(ifelse(pheno_vec_filtered == 1, "C4", "C3"))

res <- phylANOVA(tree_pruned, groups, rtt_all, nsim = 1000)
print(res)

cat("Saving results...\n")
out <- data.frame(
  gene = names(pvals),
  n_tips = ntips,
  p = pvals,
  q = qvals,
  rtt_diff = round(rtt_diff, 4),
  direction = ifelse(rtt_diff > 0, "longer_in_C4", "longer_in_C3")
)

write.csv(out, "root_to_tip_C3vsC4.csv", row.names = FALSE, quote = TRUE)
writeLines(out$gene[out$q < 0.05 & out$direction == "longer_in_C3"], "genes_C3.txt")
writeLines(out$gene[out$q < 0.05 & out$direction == "longer_in_C4"], "genes_C4.txt")

cat("\nDirection of significant genes (q < 0.05):\n")
print(table(subset(out, q < 0.05)$direction))

cat("Total genes analyzed:\n")
print(nrow(out))
cat("Significant genes (q < 0.05):\n")
print(sum(out$q < 0.05))

cat("Generating plots...\n")
dir.create("plots", showWarnings = FALSE)

# Volcano plot
volcano <- ggplot(out, aes(x = rtt_diff, y = -log10(q), color = q < 0.05)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot: Root-to-Tip C4 vs C3", x = "Mean RTT (C4 - C3)", y = "-log10(q-value)")
ggsave("plots/volcano_plot.png", volcano, width = 8, height = 6)

# Top significant genes boxplots
top_genes <- head(out[order(out$q), ], 9)$gene
box_data <- data.frame()

for (g in top_genes) {
  dna <- read.dna(file.path("gene_fastas", g), "fasta")
  rownames(dna) <- clean_names(rownames(dna))
  spp <- intersect(rownames(dna), names(pheno_vec))
  dmat <- dist.dna(dna[spp, ], "JC69")
  gtree <- njs(dmat)
  rtt <- node.depth.edgelength(gtree)[match(gtree$tip.label, gtree$tip.label)]
  names(rtt) <- gtree$tip.label
  df <- data.frame(
    gene = g,
    rtt = rtt,
    group = ifelse(names(rtt) %in% names(pheno_vec[pheno_vec == 1]), "C4", "C3")
  )
  box_data <- rbind(box_data, df)
}

p_box <- ggplot(box_data, aes(x = group, y = rtt, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "RTT by group for top genes", x = "Group", y = "Root-to-Tip Distance")
ggsave("plots/boxplots_top_genes.png", p_box, width = 10, height = 6)

# Histograms
p_hist <- ggplot(out, aes(x = p)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of p-values", x = "p-value", y = "Count")
ggsave("plots/hist_pvalues.png", p_hist, width = 6, height = 4)

q_hist <- ggplot(out, aes(x = q)) +
  geom_histogram(bins = 30, fill = "salmon", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of q-values", x = "q-value", y = "Count")
ggsave("plots/hist_qvalues.png", q_hist, width = 6, height = 4)

# Tree with phenotype annotation
colors <- ifelse(tree$tip.label %in% names(pheno_vec),
                 ifelse(pheno_vec[tree$tip.label] == 1, "firebrick", "dodgerblue"),
                 "gray80")

png("plots/tree_C3_C4.png", width = 1000, height = 1000)
plot(tree, tip.color = colors, cex = 0.6)
legend("bottomleft", legend = c("C3", "C4", "Other"),
       fill = c("dodgerblue", "firebrick", "gray80"), bty = "n")
dev.off()

cat("\nResults saved to 'root_to_tip_C3vsC4.csv' and 'plots/' folder\n")
cat("Significant gene IDs saved to 'genes_C3.txt' and 'genes_C4.txt'\n")
