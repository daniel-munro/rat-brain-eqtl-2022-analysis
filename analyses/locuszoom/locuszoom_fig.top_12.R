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


eqtls <- read_tsv("data/eqtls/top_assoc.txt", col_types = "ccciiciiccdddddddd") |>
    filter(tissue == "NAcc")

# Get samples to subset genotypes when calculating LD:
samples <- read_tsv("data/samples.txt", col_types = "cc") |>
    separate(library, c("rat_id", "old_tissue")) |>
    filter(brain_region == "NAcc")

# Top 12 eQTLs
top <- eqtls |>
    arrange(pval_beta) |>
    slice(1:12)

# Load all gene-variant pairs for them and plot.
pairs <- top |>
    mutate(gene_id = fct_inorder(gene_id)) |>
    group_by(gene_id, gene_name) |>
    summarise(
        read_tsv(str_glue("data/tensorqtl/genes/AQCT.{gene_id}.txt.gz"),
                 col_types = "-ci---d--"),
        .groups = "drop"
    ) |>
    separate(variant_id, c("chrom", "pos"), sep=":", convert = TRUE,
             remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", "") |> as.integer()) |>
    group_by(gene_id) |>
    mutate(top = pval_nominal == min(pval_nominal),
           LD = LD_with_top(variant_id, chrom, pos,
                            variant_id[top][ceiling(sum(top) / 2)],
                            samples$rat_id)) |>
    ungroup() |>
    mutate(label = str_glue("{gene_name} (chr{chrom})") |>
               fct_inorder())

egene_stats <- pairs |>
    mutate(tss = (pos - tss_distance) / 1e6) |>
    distinct(gene_id, tss) |>
    left_join(select(eqtls, gene_id, gene_name, chrom, pval_nominal_threshold),
              by = c("gene_id")) |>
    mutate(log10_threshold = -log10(pval_nominal_threshold),
           label = str_glue("{gene_name} (chr{chrom})") |>
               factor(levels = levels(pairs$label)))

ranges <- egene_stats |>
    group_by(label) |>
    summarise(pos = c(tss - 1, tss + 1), .groups = "drop") |>
    mutate(log10p = 0, LD = 0)

pairs |>
    mutate(log10p = -log10(pval_nominal),
           pos = pos / 1e6) |>
    ggplot(aes(x = pos, y = log10p, color = LD)) +
    facet_wrap(~ label, ncol = 3, scales = "free_x") +
    geom_hline(aes(yintercept = log10_threshold), data = egene_stats,
               lty = "12", size = 0.5, color = "#555555") +
    geom_vline(aes(xintercept = tss), data = egene_stats) +
    geom_point(size = 0.3) +
    expand_limits(y = 0) +
    geom_blank(data = ranges) +
    scale_x_continuous(breaks = 1:500) +
    scale_color_viridis_c(direction = -1, option = "C") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        # strip.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        panel.spacing.x = unit(10, "pt"), # to prevent tick label overlap
        legend.position = "top",
        legend.margin = margin(0, 0, -5, 0),
    ) +
    xlab("Position (Mb)") +
    ylab(expression(-log[10](P))) +
    labs(color = expression("LD "*(r^2)))

ggsave("analysis/locuszoom/top_eqtls.png", width = 5, height = 6)
