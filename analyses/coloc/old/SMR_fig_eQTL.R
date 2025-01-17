library(tidyverse)

chr_len <- c(282763074, 266435125, 177699992, 184226339, 173707219,
             147991367, 145729302, 133307652, 122095297, 112626471,
             90463843, 52716770, 114033958, 115493446, 111246239,
             90668790, 90843779, 88201929, 62275575, 56205956)
label_locs <- cumsum(c(0, chr_len[1:19])) + chr_len / 2

genes <- read_tsv("data/genes.txt", col_types = "cc----------")

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "cc---c-----------")

smr_snps <- read_tsv("coloc/SMR.tsv", col_types = "cccdcddd") |>
    filter(tissue == "PL",
           trait == "retrofat")

esnps <- read_tsv("data/tensorqtl/PQCT.cis_qtl_signif.txt.gz", col_types = "cc----ddd-") |>
    rename(gene_id = phenotype_id) |>
    filter(gene_id %in% eqtls$gene_id) |>
    left_join(
        smr_snps |>
            select(gene_id, variant_id) |>
            mutate(tested = TRUE),
        by = c("gene_id", "variant_id")
    ) |>
    replace_na(list(tested = FALSE))

df <- esnps |>
    # TMP
    slice_sample(n = 100000) |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(pval_nominal)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(gpos = pos + cumsum(c(0, chr_len))[chrom]) |>
    arrange(tested, gpos)

ggplot(df, aes(x = gpos, y = logp, color = tested)) +
    geom_point(size = 0.5) +
    expand_limits(x = c(0, sum(chr_len))) +
    expand_limits(y = c(0, max(df$logp) * 1.03)) +
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_color_manual(values = c("#aaaaaa", "black")) +
    theme_bw() +
    theme(
        legend.position = c(0.98, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_rect(color = "black", size = 0.2),
        panel.grid = element_blank(),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10](P[eQTL]))) +
    labs(color = "Tested for\ncolocalization")

ggsave("coloc/SMR_fig_eQTL.png", width = 7, height = 2)
