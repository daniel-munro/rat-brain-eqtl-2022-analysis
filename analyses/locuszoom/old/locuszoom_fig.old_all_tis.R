suppressPackageStartupMessages(library(tidyverse))

load_geno <- function(chrom, start, end, samples) {
    filename <- "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
    rng <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))
    gt <- VariantAnnotation::readGT(filename, param = rng)
    geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)[x])
    rownames(geno) <- rownames(gt)
    t(geno)[samples, ]
}

LD_with_top <- function(variant_id, chrom, pos, top_id, samples) {
    geno <- load_geno(chrom[1], min(pos) - 1, max(pos) + 1, samples)
    map_dbl(variant_id, ~ if (sd(geno[, .x]) == 0) {0} else {cor(geno[, .x], geno[, top_id]) ^ 2})
}


top_assoc <- read_tsv("data/eqtls/top_assoc.txt", col_types = "ccciiciiccdddddddd")

# Get samples to subset genotypes when calculating LD:
samples <- read_tsv("data/samples.txt", col_types = "cc") |>
    separate(library, c("rat_id", "old_tissue"))

# # Determine top ~3 genes as those with lowest maximum top-variant p-value.
# top <- eqtls |>
#     group_by(gene_id, chrom) |>
#     filter(n() == 5) |>
#     summarise(max_p = max(pval), .groups = "drop") |>
#     arrange(max_p)
# # Top per tissue
# top <- eqtls |>
#     group_by(tissue) |>
#     arrange(pval) |>
#     slice(1) |>
#     ungroup()
# Overall top 12, keeping only first instance per gene:
top <- top_assoc |>
    arrange(pval_beta) |>
    group_by(gene_id) |>
    slice(1) |>
    ungroup() |>
    arrange(pval_beta) |>
    slice(1:12)

# Load all gene-variant pairs for them and plot.

# codes <- c(Acbc = "AQCT", IL = "IQCT", LHB = "LQCT", PL = "PQCT", VoLo = "VQCT")
# pairs <- crossing(tissue = unique(eqtls$tissue),
#                   gene_id = top$gene_id[1:3]) |>
pairs <- top |>
    group_by(tissue, gene_id) |>
    summarise(
        read_tsv(str_glue("data/tensorqtl/genes/{str_sub(tissue, 1, 1)}QCT.{gene_id}.txt.gz"),
                 col_types = "-ci---d--"),
        .groups = "drop"
    ) |>
    separate(variant_id, c("chrom", "pos"), sep=":", convert = TRUE,
             remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", "") |> as.integer()) |>
    group_by(tissue, gene_id) |>
    mutate(top = pval_nominal == min(pval_nominal),
           LD = LD_with_top(variant_id, chrom, pos,
                            variant_id[top][ceiling(sum(top) / 2)],
                            samples$rat_id[samples$brain_region == unique(tissue)])) |>
    ungroup() |>
    left_join(distinct(top_assoc, gene_id, gene_name), by = "gene_id")

egene_stats <- pairs |>
    mutate(tss = (pos - tss_distance) / 1e6) |>
    distinct(tissue, gene_id, tss) |>
    left_join(select(top_assoc, tissue, gene_id, gene_name, chrom, pval_nominal_threshold),
              by = c("tissue", "gene_id")) |>
    mutate(log10_threshold = -log10(pval_nominal_threshold),
           label = str_glue("{gene_name} in {tissue} (chr{chrom})"))

ranges <- egene_stats |>
    group_by(label) |>
    summarise(pos = c(tss - 1, tss + 1), .groups = "drop") |>
    mutate(log10p = 0, LD = 0)

pairs |>
    mutate(log10p = -log10(pval_nominal),
           pos = pos / 1e6,
           label = str_glue("{gene_name} in {tissue} (chr{chrom})")) |>
    ggplot(aes(x = pos, y = log10p, color = LD)) +
    facet_wrap(~ label, ncol = 3, scales = "free_x") +
    geom_hline(aes(yintercept = log10_threshold), data = egene_stats,
               lty = "12", size = 0.5, color = "#555555") +
    geom_vline(aes(xintercept = tss), data = egene_stats) +
    geom_point(size = 0.3) +
    expand_limits(y = 0) +
    geom_blank(data = ranges) +
    scale_color_viridis_c(direction = -1, option = "C") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        # strip.background = element_rect(fill = "white"),
        strip.background = element_blank(),
    ) +
    xlab("Position (Mb)") +
    ylab(expression(-log[10]*P))

ggsave("locuszoom/top_eqtls.png", width = 8, height = 6)
