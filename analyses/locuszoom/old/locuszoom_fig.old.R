suppressPackageStartupMessages(library(tidyverse))

load_geno <- function(chrom, start, end, samples) {
    filename <- "../data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
    rng <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))
    gt <- VariantAnnotation::readGT(filename, param = rng)
    geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)[x])
    rownames(geno) <- rownames(gt)
    t(geno)[samples, ]
}

LD_with_top <- function(variant_id, chrom, pos, top_id, samples) {
    geno <- load_geno(chrom[1], min(pos) - 1, max(pos) + 1, samples)
    map_dbl(variant_id, ~ if (sd(geno[, .x]) == 0) 0 else cor(geno[, .x], geno[, top_id]) ^ 2)
}


eqtls <- read_tsv("../data/eqtls.txt", col_types = "cc-i-ddd--")
# Get samples to subset genotypes when calculating LD:
# samples <- eqtls %>%
#     group_by(tissue, region) %>%
#     summarise(
#         read_tsv(str_glue("../data/expression/{tissue}.rsem_TPM.bed.gz"),
#                  col_types = cols(`#chr` = "-", start = "-", end = "-",
#                                   gene_id = "-", .default = "d"),
#                  n_max = 1) %>%
#             pivot_longer(everything(), names_to = "sample") %>%
#             select(sample),
#         .groups = "drop"
#     )
samples <- tibble(tissue = c("Acbc", "IL", "LHB", "PL", "VoLo")) %>%
    group_by(tissue) %>%
    summarise(
        library = read.csv(str_glue("~/Dropbox (Scripps Research)/HS-RNASeq/quantitation/EnsemblGene_v2/log2+1/ensembl-gene_log2_{tissue}.txt"),
                          sep = "\t", check.names = FALSE, nrows = 1) %>%
            colnames(),
        .groups = "drop"
    ) %>%
    separate(library, c("sample", "tissue2"))

# Determine top ~3 genes as those with lowest maximum top-variant p-value.
top <- eqtls %>%
    group_by(gene_id, chrom) %>%
    filter(n() == 5) %>%
    summarise(max_p = max(pval), .groups = "drop") %>%
    arrange(max_p)

# Load all gene-variant pairs for them and plot.

codes <- c(Acbc = "AQCT", IL = "IQCT", LHB = "LQCT", PL = "PQCT", VoLo = "VQCT")
pairs <- crossing(tissue = unique(eqtls$tissue),
                  gene_id = top$gene_id[1:3]) %>%
    group_by(tissue, gene_id) %>%
    summarise(
        read_tsv(str_glue("../data/tensorqtl/genes/{codes[tissue]}.{gene_id}.txt.gz"),
                 col_types = "-ci---d--"),
        .groups = "drop"
    ) %>%
    separate(variant_id, c("chrom", "pos"), sep=":", convert = TRUE,
             remove = FALSE) %>%
    mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer()) %>%
    group_by(tissue, gene_id) %>%
    mutate(top = pval_nominal == min(pval_nominal),
           LD = LD_with_top(variant_id, chrom, pos,
                            variant_id[top][ceiling(sum(top) / 2)],
                            samples$sample[samples$tissue == unique(tissue)])) %>%
    ungroup()

egene_stats <- pairs %>%
    mutate(tss = (pos - tss_distance) / 1e6) %>%
    distinct(gene_id, tss) %>%
    left_join(select(eqtls, tissue, gene_id, pval_nominal_threshold), by = "gene_id") %>%
    mutate(log10_threshold = -log10(pval_nominal_threshold))

pairs %>%
    mutate(log10p = -log10(pval_nominal),
           LD = case_when(
               LD >= 0.99 ~ "Top",
               LD >= 0.8 ~ "[0.8 - 0.99)",
               LD >= 0.6 ~ "[0.6 - 0.8)",
               LD >= 0.4 ~ "[0.4 - 0.6)",
               LD >= 0.2 ~ "[0.2 - 0.4)",
               TRUE ~ "[0 - 0.2)"
           ) %>% fct_rev(),
           pos = pos / 1e6) %>%
    ggplot(aes(x = pos, y = log10p, color = LD)) +
    facet_grid(rows = vars(tissue), cols = vars(gene_id), scales = "free") +
    geom_hline(aes(yintercept = log10_threshold), data = egene_stats,
               lty = 2, color = "#555555") +
    geom_vline(aes(xintercept = tss), data = egene_stats) +
    geom_point(size = 0.3) +
    theme_bw() +
    xlab("Position (Mbp)") +
    ylab("-log10(p-value)")

ggsave("top_eqtls.png", width = 8, height = 6)
