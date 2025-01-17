library(tidyverse)
library(patchwork)

labels <- c("5' UTR", "3' UTR", "Missense", "Synonymous", "Splice region", "Intron",
            "Noncoding transcript", "Intergenic - upstream", "Intergenic - downstream",
            "Intergenic - other")

# Include all significant, or all tied top SNPs? The latter would exclude
# independent eQTLs, but maybe that's fine for proportions.
all_top <- tibble(tissue = c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("data/tensorqtl/{str_sub(tissue, 1, 1)}QCT.cis_qtl_signif.txt.gz"),
                 col_types = "cc----d---"),
        .groups = "drop"
    ) |>
    group_by(tissue, phenotype_id) |>
    filter(pval_nominal == min(pval_nominal)) |>
    ungroup() |>
    distinct(tissue, variant_id)

vep <- read_tsv("vep/processed_vep.txt.gz", col_types = "cc-") |>
    mutate(
        esnp_IL = variant_id %in% all_top$variant_id[all_top$tissue == "IL"],
        esnp_LHb = variant_id %in% all_top$variant_id[all_top$tissue == "LHb"],
        esnp_NAcc = variant_id %in% all_top$variant_id[all_top$tissue == "NAcc"],
        esnp_OFC = variant_id %in% all_top$variant_id[all_top$tissue == "OFC"],
        esnp_PL = variant_id %in% all_top$variant_id[all_top$tissue == "PL"]
    )

snp_counts <- vep |>
    distinct(variant_id, .keep_all = TRUE) |>
    with(c(IL = sum(esnp_IL),
           LHb = sum(esnp_LHb),
           NAcc = sum(esnp_NAcc),
           OFC = sum(esnp_OFC),
           PL = sum(esnp_PL)))

enrich <- vep |>
    group_by(Consequence) |>
    summarise(IL = sum(esnp_IL),
              LHb = sum(esnp_LHb),
              NAcc = sum(esnp_NAcc),
              OFC = sum(esnp_OFC),
              PL = sum(esnp_PL),
              n_total = n(),
              .groups = "drop") |>
    pivot_longer(IL:PL, names_to = "tissue", values_to = "n_esnp") |>
    mutate(
        #frac_esnp = n_esnp / sum(n_esnp),
        #frac_total = n_total / sum(n_total),
        frac_esnp = n_esnp / snp_counts[tissue],
        frac_total = n_total / n_distinct(vep$variant_id),
        log2_enrich = log2(frac_esnp / frac_total)
    ) |>
    mutate(Consequence = fct_relevel(Consequence, labels)) |>
    arrange(Consequence)

p1 <- enrich |>
    group_by(Consequence) |>
    summarise(mean_enr = mean(log2_enrich),
              sd_enr = sd(log2_enrich)) |>
    mutate(Consequence = fct_rev(Consequence)) |>
    ggplot(aes(x = Consequence, y = mean_enr, ymin = mean_enr - sd_enr, ymax = mean_enr + sd_enr)) +
    geom_pointrange(fatten = 1) +
    geom_hline(yintercept = 0, lty = 2) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab(NULL) +
    ylab(expression(log[2]*" fold enrichment"))

# ggsave("vep/vep_figure_enrich.png", width = 3.5, height = 2)

p2 <- enrich |>
    group_by(Consequence) |>
    summarise(mean_frac = mean(frac_esnp),
              sd_frac = sd(frac_esnp)) |>
    mutate(Consequence = fct_rev(Consequence)) |>
    ggplot(aes(x = Consequence, y = mean_frac, ymin = mean_frac - sd_frac, ymax = mean_frac + sd_frac)) +
    # geom_col(fill = "#aaaaaa") +
    geom_col(color = "black", fill = "white", size = 0.25, width = 0.5) +
    geom_linerange() +
    coord_flip() +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank()) +
    xlab(NULL) +
    ylab("Proportion\nof variants")

# ggsave("vep/vep_figure_prop.png", width = 1.5, height = 2)

p1 + p2 + plot_layout(widths = c(3, 2))
ggsave("vep/vep_figure.png", width = 5, height = 2.2)

##################################
## Effect size of each category ##
##################################

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    left_join(select(vep, variant_id, Consequence), by = "variant_id")

eqtls |>
    mutate(Consequence = Consequence |>
               fct_relevel(labels) |>
               fct_rev()) |>
    ggplot(aes(x = Consequence, y = abs(log2_aFC))) +
    geom_boxplot(outlier.size = 0.25) +
    coord_flip(ylim = c(0, 2)) +
    # coord_cartesian(xlim = c(NA, 2)) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank())
