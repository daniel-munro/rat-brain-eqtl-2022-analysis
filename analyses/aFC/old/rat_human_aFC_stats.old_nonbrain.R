# Based on rat_human_aFC.Rmd

library(tidyverse)

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciciiidddddidd") |>
    filter(rank == 1) |>
    group_by(gene_id) |>
    summarise(mean_abs_rat = mean(abs(log2_aFC)))

# # Get homologs
# # (https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/)
# human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# rat <- biomaRt::useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
# gene_map <- biomaRt::getLDS(
#     attributes = c("ensembl_gene_id"),
#     filters = "ensembl_gene_id",
#     values = unique(eqtls$gene_id),
#     mart = rat,
#     attributesL = c("ensembl_gene_id_version"),
#     martL = human,
#     uniqueRows = TRUE
# ) |>
#     rename(gene_id_rat = Gene.stable.ID,
#            gene_id_human = Gene.stable.ID.version)
# write_tsv(gene_map, "aFC/gene_map.txt")
gene_map <- read_tsv("aFC/gene_map.txt", col_types = "cc")

gtex <- list.files("data/gtex/GTEx_Analysis_v8_eQTL", full.names = TRUE) |>
    read_tsv(id = "file",
             col_types = cols(gene_id = "c", qval = "d", pval_beta = "-",
                              log2_aFC = "d", .default = "-")) |>
    filter(qval <= 0.05) |>  # GTEx site says to do <=
    mutate(group = if_else(str_detect(file, "Brain"), "GTEx_brain", "GTEx_other")) |>
    group_by(group, gene_id) |>
    summarise(mean_abs_human = mean(abs(log2_aFC)),
              .groups = "drop")

# All homolog pairs with eQTLs in rat and human brain, same with human non-brain:
d <- gene_map |>
    inner_join(eqtls, by = c("gene_id_rat" = "gene_id")) |>
    inner_join(gtex, by = c("gene_id_human" = "gene_id"))

stats <- d |>
    group_by(group) |>
    summarise(
        pairs = n(),
        pearson = cor(mean_abs_rat, mean_abs_human),
        pearson_p = cor.test(mean_abs_rat, mean_abs_human)$p.value,
        spearman = cor(mean_abs_rat, mean_abs_human, method = "spearman"),
        spearman_p = cor.test(mean_abs_rat, mean_abs_human, method = "spearman")$p.value,
        .groups = "drop"
    )

deming <- d |>
    mutate(group = fct_recode(group,
                              "Human brain tissues" = "GTEx_brain",
                              "Human non-brain tissues" = "GTEx_other")) |>
    group_by(group) |>
    summarise({
        coef <- deming::deming(mean_abs_human ~ mean_abs_rat)$coefficients
        tibble(intercept = coef[1],
               slope = coef[2])
    })

d |>
    mutate(group = fct_recode(group,
                              "Human brain tissues" = "GTEx_brain",
                              "Human non-brain tissues" = "GTEx_other")) |>
    ggplot(aes(x = mean_abs_rat, y = mean_abs_human)) +
    facet_wrap(~ group) +
    geom_point(size = 1, alpha = 0.5) +
    geom_abline(aes(intercept = intercept, slope = slope), data = deming,
                color = "#8888cc") +
    coord_fixed() +
    xlab(expression("Mean |"*log[2]*"aFC| in Rat")) +
    ylab(expression("Mean |"*log[2]*"aFC| in Human")) +
    theme_minimal()

ggsave("aFC/rat_human_aFC_fig.png", width = 5.2, height = 3, bg = "white")
