library(tidyverse)

# Look at gene strand to orient TSS distance.
genes <- read_tsv("data/genes.txt", col_types = "c----c------")

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    left_join(genes, by = "gene_id") |>
    # mutate(tss_distance = pos - tss)
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos) / 1e6)

eqtls |>
    ggplot(aes(x = tss_distance)) +
    geom_histogram(bins = 50) +
    # geom_density(alpha = 0.5) +
    # geom_vline(aes(xintercept = mean(eqtls$tss_distance)), alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0), limits = c(-1.001, 1.001)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("Distance from TSS (Mb)") +
    ylab(NULL) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())

ggsave("analysis/stats/TSS_distance.png", width = 3, height = 2.5)

#####################
## GTEx comparison ##
#####################
# Compare top eQTL per gene, not all independent eQTLs, since non-primary eQTLs
# are probably further away and there are many more of them in GTEx.

gtex_top <- tibble(
    file = list.files("data/gtex/GTEx_Analysis_v8_eQTL",
    full.names = TRUE
)) |>
    group_by(file) |>
    summarise(
        read_tsv(
            file,
            col_types = cols(
                gene_id = "c", qval = "d", pval_beta = "-",
                log2_aFC = "d", tss_distance = "i", .default = "-"
            )
        ),
        .groups = "drop"
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
    # geom_histogram(bins = 50)
    geom_density(adjust = 0.5, alpha = 0.5)

tss |>
    mutate(group = fct_rev(group),
           abs_tss_distance = abs_tss_distance / 1e6) |>
    ggplot(aes(x = abs_tss_distance, color = group)) +
    geom_density(adjust = 0.5, size = 0.8, outline.type = "full", show.legend = FALSE) +
    scale_x_continuous(expand = c(0.02, 0)) +
    # scale_y_continuous(expand = c(0.02, 0), breaks = c(0, 50, 100, 150) / max_g,
    #                    sec.axis = sec_axis(~ ., breaks = c(0, 1, 2) / max_r)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    # from (scales::hue_pal())(3)[c(1, 3)]:
    scale_color_manual(values = c("#F8766DFF", "#619CFFFF")) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        # panel.border = element_blank(),
        # axis.line.x = element_line(),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
    ) +
    xlab("TSS distance (Mb)") +
    ylab("Density") # Line up with LD decay panel

ggsave("analysis/stats/TSS_distance_compare.png", width = 2.5, height = 2.5, dpi = 300)

##############################
## Stratify by VEP category ##
##############################

labels <- c("5' UTR", "3' UTR", "Missense", "Synonymous", "Splice region", "Intron",
            "Noncoding transcript", "Intergenic - upstream", "Intergenic - downstream",
            "Intergenic - other")

vep <- read_tsv("analysis/vep/processed_vep.txt.gz", col_types = "ccc") |>
    mutate(Consequence = factor(Consequence, levels = labels))

eqtls_vep <- eqtls |>
    inner_join(vep, by = "variant_id") |>
    mutate(
        anno_gene = case_when(
            is.na(Gene) ~ "None",
            gene_id == Gene ~ "eGene",
            gene_id != Gene ~ "Nearby gene",
            TRUE ~ "Error"
        ), #|>
            # fct_relevel("Same"),
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
    # ggtitle("eQTLs colored by whether eGene matches annotation gene")

ggsave("analysis/stats/TSS_distance.VEP.png", width = 6, height = 6)

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
# ggtitle("eQTLs colored by whether eGene matches annotation gene")

ggsave("analysis/stats/TSS_distance.VEP2.png", width = 6, height = 7)

# Alternative: flip coloring and faceting

eqtls_vep |>
    ggplot(aes(x = tss_distance, fill = Consequence)) +
    facet_wrap(~ anno_gene, scales = "free_y") +
    geom_histogram(bins = 80) +
    xlab("Distance from TSS (Mb)") +
    ylab(NULL) +
    theme_minimal()
# ggtitle("eQTLs colored by whether eGene matches annotation gene")

#########
## Old ##
#########

## Originally I included tied SNP start and end distributions:
# codes <- c(Acbc = "AQCT", IL = "IQCT", LHB = "LQCT", PL = "PQCT", VoLo = "VQCT")
# 
# eqtls <- tibble(tissue = names(codes)) |>
#     group_by(tissue) |>
#     summarise(
#         read_tsv(str_glue("data/tensorqtl/{codes[tissue]}.cis_qtl.txt.gz"),
#         col_types = cols(phenotype_id = "c", variant_id = "c",
#                          tss_distance = "i", .default = "-")),
#         .groups = "drop"
#     ) |>
#     rename(gene_id = phenotype_id)
# eqtls2 <- read_tsv("data/brain.top_assoc_per_gene.txt.gz",
#                    col_types = "cc-cc--d-") |>
#     filter(qval < 0.05) |>
#     left_join(eqtls, by = c("tissue", "gene_id")) |>
#     separate(first_top_variant, into = c("first_chrom", "first_pos"),
#              convert = TRUE) |>
#     separate(last_top_variant, into = c("last_chrom", "last_pos"),
#              convert = TRUE) |>
#     mutate(last_tss_distance = tss_distance + (last_pos - first_pos))
# 
# tmp <- bind_rows(
#     tibble(tss_distance = eqtls2$tss_distance / 1e6) |>
#         mutate(`Top SNP` = "First"),
#     tibble(tss_distance = eqtls2$last_tss_distance / 1e6) |>
#         mutate(`Top SNP` = "Last")
# ) |>
#     filter(!is.na(tss_distance))
# tmp_mean <- tmp |>
#     group_by(`Top SNP`) |>
#     summarise(tss_distance = mean(tss_distance),
#               .groups = "drop")
# tmp |>
#     ggplot(aes(x = tss_distance, fill = `Top SNP`)) +
#     # geom_histogram(bins = 50) +
#     geom_density(alpha = 0.5) +
#     geom_vline(aes(xintercept = tss_distance, color = `Top SNP`),
#                data = tmp_mean, alpha = 0.5) +
#     xlab("Distance from TSS (Mbp)") +
#     ylab(NULL) +
#     theme_minimal() +
#     theme(panel.grid.major.y = element_blank(),
#           panel.grid.minor.y = element_blank(),
#           axis.text.y = element_blank())
# 
# ggsave("stats/TSS_distance.png", width = 4, height = 2.5)

# ## 2D density plot?
# 
# eqtls3 |>
#     filter(!is.na(last_tss_distance)) |>
#     ggplot(aes(x = tss_distance, y = last_tss_distance)) +
#     geom_density_2d_filled()

####################################
## Distribution of top SNP ranges ##
####################################
# Maybe misleading because sometimes one p-value is slightly lower than many others.

