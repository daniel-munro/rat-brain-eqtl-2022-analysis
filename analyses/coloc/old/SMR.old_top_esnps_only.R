library(tidyverse)

set.seed(2593579)

T_SMR <- function(z_eQTL, z_GWAS) {
    (z_eQTL^2 * z_GWAS^2) / (z_eQTL^2 + z_GWAS^2)
}

p_SMR <- function(T_SMR) {
    pchisq(T_SMR, df = 1, lower.tail = FALSE)
}

#' Given a tissue name and vector of eGenes, return tibble of all tied top SNPs per gene
top_eSNPs <- function(tissue, genes) {
    signif <- str_glue("data/tensorqtl/{str_sub(tissue, 1, 1)}QCT.cis_qtl_signif.txt.gz") |>
        read_tsv(col_types = "cc----ddd-") |>
        rename(gene_id = phenotype_id) |>
        filter(gene_id %in% genes) |>
        group_by(gene_id) |>
        mutate(min_pval = min(pval_nominal)) |>
        ungroup() |>
        ## Use tolerance for floating point discrepancies etc.
        filter(pval_nominal < min_pval * 1.000001) |>
        select(-min_pval)
}

gwas_dir <- "data/coloc/adiposity_GWAS/"
gwas <- tibble(file = list.files(gwas_dir)) |>
    group_by(file) |>
    summarise(
        read_tsv(str_c(gwas_dir, file), col_types = "-c-----dd------",
                 col_names = c("variant_id", "slope", "slope_se")),
        .groups = "drop"
    ) |>
    mutate(trait = str_match(file, "allChr_physiological_(.+)\\.assoc\\.txt")[, 2], .before = 1) |>
    select(-file) |>
    mutate(z_GWAS = slope / slope_se) |>
    select(-slope, -slope_se)

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "cc---c-----------")

df <- eqtls |>
    ## Get tied top eSNPs per eGene
    distinct(tissue, gene_id) |>
    group_by(tissue) |>
    summarise(top_eSNPs(unique(tissue), gene_id),
              .groups = "drop") |>
    mutate(z_eQTL = slope / slope_se) |>
    select(-pval_nominal, -slope, -slope_se) |>
    ## Subset to tested GWAS SNPs and choose one per tissue/gene at random
    filter(variant_id %in% gwas$variant_id) |>
    group_by(tissue, gene_id) |>
    slice_sample(n = 1) |>
    ungroup() |>
    ## Join with GWAS scores and calculate SMR
    inner_join(gwas, by = "variant_id") |>
    mutate(T_SMR = T_SMR(z_eQTL, z_GWAS),
           p_SMR = p_SMR(T_SMR))

write_tsv(df, "coloc/SMR.tsv")

sig <- df |>
    group_by(tissue, trait) |>
    filter(p_SMR < 0.05 / n()) |>
    ungroup()

write_tsv(sig, "coloc/SMR_sig.tsv")
