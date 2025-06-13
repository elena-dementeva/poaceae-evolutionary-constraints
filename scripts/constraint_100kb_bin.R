library(data.table)
library(tidyverse)
library(broom)
library(qvalue)
library(ggrepel)
library(RColorBrewer)

cat("Loading input data...\n")

bindata_raw <- fread("merged_bindata.100000.tsv", sep = "\t", header = FALSE)

col_base <- c(
  "chr","start","end",
  "totalBases","conservedBases",
  "cds.distinct.sum","exon.distinct.sum",
  "gene.count","intergenic.length","GC",
  "TE.count","ISBP.count",
  "totalConservedBases","meanConsScore",
  "gene_id","gene_name","description","GO"
)

if (ncol(bindata_raw) == length(col_base)) {
  setnames(bindata_raw, col_base)
  bindata_raw$medianConsScore <- bindata_raw$meanConsScore
} else if (ncol(bindata_raw) == length(col_base) + 1) {
  setnames(bindata_raw, c(col_base, "medianConsScore"))
} else {
  stop("Unexpected number of columns in merged_bindata.100000.tsv")
}

bindata <- bindata_raw %>%
  mutate(across(
    c(start:end, totalBases:ISBP.count,
      totalConservedBases, meanConsScore, medianConsScore),
    as.numeric
  )) %>%
  filter(
    str_detect(chr, "A$"),
    (end - start) == 100000,
    !is.na(conservedBases)
  )

cat("Fitting linear model (medianConsScore predictor)…\n")
model_data <- bindata %>%
  drop_na(totalConservedBases, totalBases, cds.distinct.sum, gene.count,
          intergenic.length, GC, TE.count, ISBP.count, medianConsScore)

fit <- lm(sqrt(totalConservedBases) ~ 
            sqrt(totalBases) +
            sqrt(cds.distinct.sum + 1) +
            sqrt(gene.count + 1) +
            sqrt(intergenic.length + 1) +
            GC +
            sqrt(TE.count + 1) +
            sqrt(ISBP.count + 1) +
            sqrt(medianConsScore + 1),
          data = model_data)

aug <- augment(fit, model_data) %>%
  mutate(pval = 2 * pnorm(abs(.std.resid), lower.tail = FALSE))

aug$qval <- qvalue(aug$pval)$qvalues
aug$sig  <- aug$qval <= 0.05

aug_full <- bindata %>%
  left_join(
    as.data.frame(aug) %>%
      dplyr::select(chr, start, end, .fitted, .resid, .std.resid, pval, qval, sig),
    by = c("chr", "start", "end")
  )

cat("Preparing Manhattan coordinates…\n")
gap_size <- 2e6
chrom_sizes <- aug_full %>%
  filter(!is.na(qval)) %>%
  group_by(chr) %>%
  summarise(chr_len = max(end), .groups = "drop") %>%
  arrange(chr) %>%
  mutate(chr_start = lag(cumsum(chr_len + gap_size), default = 0))

aug_annot <- aug_full %>%
  left_join(chrom_sizes, by = "chr") %>%
  mutate(
    pos_cum    = start + chr_start,
    gene_label = ifelse(is.na(gene_name) | gene_name == "NA", gene_id, gene_name),
    is_auto    = gene_label == gene_id
  )

cat("Assigning numeric labels…\n")
auto_ids <- aug_annot %>%
  filter(sig, is_auto) %>%
  distinct(gene_id) %>%
  mutate(num_label = as.character(row_number()))

aug_annot <- aug_annot %>%
  left_join(auto_ids, by = "gene_id") %>%
  mutate(gene_plot_label = if_else(is_auto, num_label, gene_label)) %>%
  as_tibble()

cat("Building legend for numeric labels…\n")
label_legend <- aug_annot %>%
  filter(!is.na(num_label)) %>%
  distinct(num_label, gene_id, description, GO) %>%
  mutate(
    desc_clean = ifelse(description %in% c(NA, "", "NA"), "", description),
    go_clean   = ifelse(GO         %in% c(NA, "", "NA"), "", GO),
    legend_line = case_when(
      desc_clean != "" & go_clean != "" ~ paste0(num_label, " - ", gene_id,
                                                 " — ", desc_clean, " [", go_clean, "]"),
      desc_clean != "" & go_clean == "" ~ paste0(num_label, " - ", gene_id,
                                                 " — ", desc_clean),
      desc_clean == "" & go_clean != "" ~ paste0(num_label, " - ", gene_id,
                                                 " [", go_clean, "]"),
      TRUE                               ~ paste0(num_label, " - ", gene_id)
    )
  ) %>%
  arrange(as.numeric(num_label))

writeLines(label_legend$legend_line, "significant_gene_legend.txt")
cat("Legend saved.\n")

cat("Selecting top windows per chromosome…\n")
sig_top <- aug_annot %>%
  filter(sig) %>%
  group_by(chr) %>%
  slice_max(order_by = medianConsScore, n = 10) %>%
  ungroup()

color_vec  <- brewer.pal(min(12, length(unique(aug_annot$chr))), "Set3")
names(color_vec) <- unique(aug_annot$chr)

cat("Plotting –log10(q)…\n")

plot_title <- "100 kb Constraint Manhattan Plot (A subgenome)"
manh_q <- ggplot() +
  geom_point(data = aug_annot,
             aes(pos_cum, -log10(qval), colour = chr),
             size = 0.6, alpha = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = sig_top,
                  aes(pos_cum, -log10(qval), label = gene_plot_label),
                  size = 2.4, box.padding = 0.4, max.overlaps = 100) +
  scale_colour_manual(values = color_vec, guide = "none") +
  scale_x_continuous(breaks = chrom_sizes$chr_start + chrom_sizes$chr_len/2,
                     labels = chrom_sizes$chr) +
  labs(x = "Chromosome", y = expression(-log[10](qvalue)), title = plot_title) +
  theme_bw()
ggsave("100kb_constraint_manhattan_q.png", manh_q, width = 12, height = 5, dpi = 300)

cat("Plotting medianConsScore…\n")

sig_top_med <- aug_annot %>%
  filter(sig) %>%
  group_by(chr) %>%
  slice_max(order_by = medianConsScore, n = 10) %>%
  ungroup()

manh_med <- ggplot() +
  geom_point(data = aug_annot,
             aes(pos_cum, medianConsScore, colour = chr),
             size = 0.6, alpha = 0.7) +
  geom_text_repel(data = sig_top_med,
                  aes(pos_cum, medianConsScore, label = gene_plot_label),
                  size = 2.4, box.padding = 0.4, max.overlaps = 100) +
  scale_colour_manual(values = color_vec, guide = "none") +
  scale_x_continuous(breaks = chrom_sizes$chr_start + chrom_sizes$chr_len / 2,
                     labels = chrom_sizes$chr) +
  labs(x = "Chromosome", y = "medianConsScore",
       title = "100 kb Constraint Manhattan Plot (A subgenome)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("100kb_constraint_manhattan_medianScore.png", manh_med, width = 12, height = 5, dpi = 300)


