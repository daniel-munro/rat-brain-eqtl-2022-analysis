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

# Top 6 eQTLs
top <- eqtls |>
    arrange(pval_beta) |>
    slice(1:6)

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
        panel.spacing.y = unit(90, "pt"), # to make space for effect plots
        legend.position = "top",
        legend.margin = margin(0, 0, -5, 0),
    ) +
    xlab(NULL) +
    ylab(expression(-log[10](P))) +
    labs(color = expression("LD "*(r^2)))

ggsave("analysis/locuszoom/top_eqtls_locus.png", width = 5, height = 4.5)

##################
## Effect plots ##
##################

library(patchwork)

load_snp <- function(chrom, pos) {
    filename <- "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
    rng <- GenomicRanges::GRanges(chrom, IRanges::IRanges(pos, pos))
    gt <- VariantAnnotation::readGT(filename, param = rng)
    apply(gt, 2, function(x) c("0|0" = 0L, "0|1" = 1L, "1|0" = 1L, "1|1" = 2L)[x])
}

plot_eqtl <- function(df) {
    df |>
        arrange(geno) |>
        rowwise() |>
        mutate(
            geno = c(str_glue("{ref}/{ref}"),
                     str_glue("{ref}/{alt}"),
                     str_glue("{alt}/{alt}"))[geno + 1] |>
                fct_inorder()
        ) |>
        ungroup() |>
        ggplot(aes(x = geno, y = expr, fill = geno)) +
        # facet_wrap(~ variant_id, ncol = 3, scales = "free_x") +
        geom_violin(show.legend = FALSE) +
        geom_boxplot(fill = "white", width = 0.1, outlier.size = 0.5) +
        scale_fill_manual(values = c("#ebfbe5", "#95e879", "#4ac021")) +
        xlab(NULL) +
        ylab(NULL) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            # strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 10),
        ) +
        ggtitle(unique(df$variant_id))
}

expr <- read_tsv("data/expression/ensembl-gene_inv-quant_NAcc.bed.gz",
                 col_types = cols(`#chr` = "-", start = "-", end = "-",
                                  gene_id = "c", .default = "d")) |>
    pivot_longer(-gene_id, names_to = "rat_id", values_to = "expr")

effect <- top |>
    group_by(gene_id, gene_name, variant_id, chrom, pos, ref, alt) |>
    summarise(
        load_snp(chrom, pos) |>
            enframe(name = "rat_id", value = "geno") |>
            filter(rat_id %in% samples$rat_id),
        .groups = "drop"
    ) |>
    left_join(expr, by = c("gene_id", "rat_id"))

p <- lapply(top$gene_id, function (x) {
    effect |>
        filter(gene_id == x) |>
        plot_eqtl()
})
(p[[1]] + p[[2]] + p[[3]]) / plot_spacer() / (p[[4]] + p[[5]] + p[[6]]) +
    plot_layout(heights = c(2, 3, 2)) &
    theme(
        plot.background = element_blank(),
        plot.margin = margin(r = -4),
    )
ggsave("analysis/locuszoom/top_eqtls_effect.png", width = 6, height = 5)
