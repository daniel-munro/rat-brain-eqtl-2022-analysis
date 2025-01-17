library(tidyverse)

heatmap_order <- function(fac1, fac2, value) {
    df <- tibble(fac1, fac2, value) |>
        pivot_wider(fac1, names_from = fac2, values_from = value)
    mat <- df |> select(-fac1) |> as.matrix()
    rownames(mat) <- df$fac1
    clust <- hclust(dist(mat))
    factor(fac1, levels = clust$labels[clust$order])
}

tissues <- c(IL = "Infralimbic cortex",
             LHb = "Lateral habenula",
             NAcc = "Nucleus accumbens core",
             OFC = "Orbitofrontal cortex",
             PL = "Prelimbic cortex")

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciciiidddddidd") |>
    filter(rank == 1) |>
    mutate(abs_aFC_rat = abs(log2_aFC),
           tissue = tissues[tissue]) |>
    select(tissue, gene_id, abs_aFC_rat)

gene_map <- read_tsv("aFC/gene_map.txt", col_types = "cc")

gtex <- list.files("data/gtex/GTEx_Analysis_v8_eQTL", pattern = "*Brain_*",
                   full.names = TRUE) |>
    read_tsv(id = "file",
             col_types = cols(gene_id = "c", qval = "d", pval_beta = "-",
                              log2_aFC = "d", .default = "-")) |>
    filter(qval <= 0.05) |>  # GTEx site says to do <=
    mutate(tissue = str_match(file, "Brain_(.+)\\.v8")[, 2] |>
               str_replace_all("_", " "),
           abs_aFC_human = abs(log2_aFC)) |>
    select(tissue, gene_id, abs_aFC_human)

# # All homolog pairs with eQTLs in rat and human brain:
# d <- gene_map |>
#     inner_join(eqtls, by = c("gene_id_rat" = "gene_id")) |>
#     inner_join(gtex, by = c("gene_id_human" = "gene_id"))

pairs <- crossing(rat_tissue = unique(eqtls$tissue),
                  human_tissue = unique(gtex$tissue)) |>
    group_by(rat_tissue, human_tissue) |>
    summarise({
        d <- gene_map |>
            inner_join(filter(eqtls, tissue == rat_tissue),
                       by = c("gene_id_rat" = "gene_id")) |>
            inner_join(filter(gtex, tissue == human_tissue),
                       by = c("gene_id_human" = "gene_id"))
        tibble(eGene_overlap = nrow(d),
               aFC_corr_p = with(d, cor(abs_aFC_rat, abs_aFC_human)),
               aFC_corr_s = with(d, cor(abs_aFC_rat, abs_aFC_human, method = "spearman")))
    }, .groups = "drop")

pairs |>
    ggplot(aes(x = human_tissue, y = rat_tissue, fill = eGene_overlap)) +
    geom_tile() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))

pairs |>
    mutate(human_tissue = heatmap_order(human_tissue, rat_tissue, aFC_corr_p),
           rat_tissue = heatmap_order(rat_tissue, human_tissue, aFC_corr_p)) |>
    ggplot(aes(x = human_tissue, y = rat_tissue, fill = aFC_corr_p)) +
    geom_tile() +
    scale_fill_viridis_c(option = "E") +
    # scale_fill_gradient(low = "white", high = "black") +
    expand_limits(fill = c(0, 1)) +
    xlab("Human brain tissue") +
    ylab("Rat brain tissue") +
    labs(fill = expression("|"*log[2]*"aFC| corr.")) +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))

ggsave("aFC/rat_human_heatmap.png", width = 5.5, height = 3.6)

# pairs |>
#     ggplot(aes(x = human_tissue, y = rat_tissue, fill = aFC_corr_s)) +
#     geom_tile() +
#     theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
# 
# ggsave("aFC/rat_human_heatmap2.png", width = 4, height = 4)
