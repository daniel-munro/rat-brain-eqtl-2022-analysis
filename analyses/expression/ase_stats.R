library(tidyverse)

###############################################
## What % of expressed genes had ASE counts? ##
###############################################

tissue <- c("IL", "LHb", "NAcc", "OFC", "PL")

# BED files have all annotated genes, so subset to expressed genes
genes <- tibble(tissue) |>
    reframe(
        read_tsv(
            str_glue("data/expression/ensembl-gene_log2_{tissue}.bed.gz"),
            col_types = cols(gene_id = 'c', .default = '-')
        ),
        .by = tissue
    )

ase <- tibble(tissue) |>
    reframe(
        read_tsv(
            str_glue("data/phaser/{tissue}.expr_matrix.gw_phased.bed.gz"),
            col_types = cols(`#contig` = '-', start = '-', stop = '-', name = 'c', .default = 'c')
        ) |>
            rename(gene_id = name) |>
            pivot_longer(-gene_id, names_to = 'sample', values_to = 'expr'),
        .by = tissue
    ) |>
    semi_join(genes, by = c("tissue", "gene_id")) |>
    separate_wider_delim(expr, "|", names = c("count1", "count2")) |>
    mutate(count1 = as.integer(count1),
           count2 = as.integer(count2))

expr_genes_w_ase <- ase |>
    filter(count1 + count2 > 0) |>
    distinct(tissue, gene_id) |>
    count(tissue) |>
    left_join(count(genes, tissue, name = "n_expr"), by = "tissue") |>
    mutate(frac_expr_w_ase = n / n_expr)
summary(expr_genes_w_ase)

############################################
## What % of cis-eQTLs had ASE-based aFC? ##
############################################

afc <- read_tsv("data/eqtls/ASE_aFC.txt", col_types = "cccd")

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    left_join(afc, by = c("tissue", "gene_id", "variant_id"))

afc_counts <- eqtls |>
    summarise(n_eQTLs = n(),
              n_eQTLs_with_ASE = sum(!is.na(log2_aFC_ASE)),
              .by = tissue) |>
    mutate(eqtls_w_ase = n_eQTLs_with_ASE / n_eQTLs)
summary(afc_counts)
