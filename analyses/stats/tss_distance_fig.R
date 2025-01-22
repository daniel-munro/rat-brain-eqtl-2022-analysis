library(tidyverse)

# Look at gene strand to orient TSS distance.
genes <- read_tsv("data/genes.txt", col_types = "c----c------")

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    left_join(genes, by = "gene_id") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos) / 1e6)

eqtls |>
    ggplot(aes(x = tss_distance)) +
    geom_histogram(bins = 50) +
    scale_x_continuous(expand = c(0, 0), limits = c(-1.001, 1.001)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("Distance from TSS (Mb)") +
    ylab(NULL) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())

ggsave("analyses/stats/TSS_distance.png", width = 3, height = 2.5)

#####################
## GTEx comparison ##
#####################
# Compare top eQTL per gene, not all independent eQTLs, since non-primary eQTLs
# are probably further away and there are many more of them in GTEx.

gtex_top <- tibble(
    file = list.files("data/gtex/GTEx_Analysis_v8_eQTL",
    full.names = TRUE
)) |>
    reframe(
        read_tsv(
            file,
            col_types = cols(
                gene_id = "c", qval = "d", pval_beta = "-",
                log2_aFC = "d", tss_distance = "i", .default = "-"
            )
        ),
        .by = file
    ) |>
    filter(qval <= 0.05) |> # GTEx site says to do <=
    mutate(
        tissue = str_match(file, "eQTL/(.+)\\.v8")[, 2],
        group = if_else(str_detect(tissue, "Brain"), "GTEx brain", "GTEx other"),
        abs_tss_distance = abs(tss_distance)
    )

egenes <- eqtls |>
    filter(rank == 1)

tss <- bind_rows(
    gtex_top |>
        filter(group == "GTEx brain") |>
        select(tissue, group, abs_tss_distance),
    egenes |>
        mutate(group = "Rat brain",
               abs_tss_distance = abs(pos - tss)) |>
        select(tissue, group, abs_tss_distance)
)

tss |>
    group_by(group) |>
    skimr::skim(abs_tss_distance)

tss |>
    ggplot(aes(x = group, y = abs_tss_distance)) +
    geom_violin() +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA)

tss |>
    ggplot(aes(x = abs_tss_distance, fill = group)) +
    geom_density(adjust = 0.5, alpha = 0.5)

tss |>
    mutate(group = fct_rev(group),
           abs_tss_distance = abs_tss_distance / 1e6) |>
    ggplot(aes(x = abs_tss_distance, color = group)) +
    geom_density(adjust = 0.5, size = 0.8, outline.type = "full", show.legend = FALSE) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    scale_color_manual(values = c("#F8766DFF", "#619CFFFF")) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
    ) +
    xlab("TSS distance (Mb)") +
    ylab("Density") # Line up with LD decay panel

ggsave("analyses/stats/TSS_distance_compare.png", width = 2.5, height = 2.5, dpi = 300)

##############################
## Stratify by VEP category ##
##############################

labels <- c("5' UTR", "3' UTR", "Missense", "Synonymous", "Splice region", "Intron",
            "Noncoding transcript", "Intergenic - upstream", "Intergenic - downstream",
            "Intergenic - other")

vep <- read_tsv("data/vep/processed_vep.txt.gz", col_types = "ccc") |>
    mutate(Consequence = factor(Consequence, levels = labels))

eqtls_vep <- eqtls |>
    inner_join(vep, by = "variant_id", relationship = "many-to-many") |>
    mutate(
        anno_gene = case_when(
            is.na(Gene) ~ "None",
            gene_id == Gene ~ "eGene",
            gene_id != Gene ~ "Nearby gene",
            TRUE ~ "Error"
        ),
        Consequence2 = if_else(is.na(Gene) | gene_id == Gene,
                               as.character(Consequence),
                               "Nearby gene") |>
            factor(levels = c(labels, "Nearby gene"))
    )

eqtls_vep |>
    ggplot(aes(x = tss_distance, fill = anno_gene)) +
    facet_wrap(~ Consequence, scales = "free_y", ncol = 2) +
    geom_histogram(bins = 80, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = c("#1957bd", "#a5c3f3", "#888888")) +
    xlab("Distance from TSS (Mb)") +
    ylab(NULL) +
    labs(fill = "Annotation gene") +
    theme_minimal() +
    theme(legend.position = "top")

ggsave("analyses/stats/TSS_distance.VEP.png", width = 6, height = 6)

# Alternative: Group non-eGene gene annotations

eqtls_vep |>
    ggplot(aes(x = tss_distance, fill = anno_gene)) +
    facet_wrap(~ Consequence2, scales = "free_y", ncol = 2) +
    geom_histogram(bins = 80, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = c("#1957bd", "#a5c3f3", "#888888")) +
    xlab("Distance from TSS (Mb)") +
    ylab(NULL) +
    labs(fill = "Annotation gene") +
    theme_minimal() +
    theme(legend.position = "top")

ggsave("analyses/stats/TSS_distance.VEP2.png", width = 6, height = 7)

# Alternative: flip coloring and faceting

eqtls_vep |>
    ggplot(aes(x = tss_distance, fill = Consequence)) +
    facet_wrap(~ anno_gene, scales = "free_y") +
    geom_histogram(bins = 80) +
    xlab("Distance from TSS (Mb)") +
    ylab(NULL) +
    theme_minimal()
