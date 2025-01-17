library(tidyverse)

# NOTE: I'm now using a combined e/sQTL version generated in stats/eqtl_stats_fig.R.

sqtls <- read_tsv("data/splice/sqtls_indep.txt", col_types = "cciccicdddddcii")

sgenes <- sqtls |>
    distinct(tissue, group_id) |>
    rename(gene_id = group_id)

sgenes |>
    mutate(tissue = fct_infreq(tissue) |> fct_rev()) |> # to match upset plot order
    ggplot(aes(x = tissue)) +
    geom_bar() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(table(sgenes$tissue)) * 1.02)) +
    xlab(NULL) +
    # ylab("sGenes (5% FDR)") +
    ylab("sGenes") +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
          panel.grid = element_blank())

# ggsave("splice/sGene_count.png", width = 1.5, height = 1.5)
ggsave("splice/sGene_count.png", width = 1.3, height = 2.5)

tmp <- sgenes |>
    group_by(gene_id) |>
    summarise(n_tissues = n()) |>
    ungroup() |>
    count(n_tissues)
ggplot(tmp, aes(x = n_tissues, y = n, fill = n_tissues)) +
    geom_col(show.legend = FALSE) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(tmp$n) * 1.02)) +
    scale_fill_viridis_c() +
    xlab("No. tissues") +
    # ylab("sGenes (5% FDR)") +
    ylab("sGenes") +
    theme_bw() +
    theme(panel.grid = element_blank())

# ggsave("splice/sGene_count_by_n_tissues.png", width = 1.5, height = 1.5)
ggsave("splice/sGene_count_by_n_tissues.png", width = 1.3, height = 2.5)
