library(tidyverse)

chr_len <- c(282763074, 266435125, 177699992, 184226339, 173707219,
             147991367, 145729302, 133307652, 122095297, 112626471,
             90463843, 52716770, 114033958, 115493446, 111246239,
             90668790, 90843779, 88201929, 62275575, 56205956)
label_locs <- cumsum(c(0, chr_len[1:19])) + chr_len / 2

genes <- read_tsv("data/genes.txt", col_types = "cc----------")

smr <- read_tsv("coloc/SMR.tsv", col_types = "cccdcddd") |>
    left_join(genes, by = "gene_id")
sig <- smr |>
    filter(p_SMR < 0.05 / nrow(smr))

df <- sig |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(p_SMR)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(gpos = pos + cumsum(c(0, chr_len))[chrom],
           label = str_glue("{gene_name},\n{trait}")) |>
    arrange(gpos)

ggplot(df, aes(x = gpos, y = logp, color = tissue, label = label)) +
    geom_point() +
    geom_text(hjust = 0, show.legend = FALSE) +
    expand_limits(x = c(0, sum(chr_len))) +
    scale_color_manual(values = c("#377eb8", "#4daf4a", "#e41a1c", "#ff7f00", "#984ea3")) +
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0)) +
    theme_bw() +
    theme(
        legend.position = c(0.98, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_rect(color = "black", size = 0.2),
        panel.grid = element_blank(),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10](P[SMR]))) +
    labs(color = "Tissue")

ggsave("coloc/SMR_fig.svg", width = 7, height = 3)
