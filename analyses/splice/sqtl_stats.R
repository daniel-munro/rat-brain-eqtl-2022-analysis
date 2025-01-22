library(tidyverse)

# top_assoc <- read_tsv("data/splice/top_assoc_splice.txt", col_types = "ccicciidddddcidd")

sqtls <- read_tsv("data/splice/sqtls_indep.txt", col_types = "cciccicdddddcii")

sgenes <- sqtls |>
    distinct(tissue, group_id) |>
    rename(gene_id = group_id)

egenes <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    distinct(tissue, gene_id)

# Total independent sQTLs
nrow(sqtls)

# How many 2nd, 3rd, 4th eQTLs per gene
count(sqtls, rank)

# Total genes with sQTL in any tissue
n_distinct(sgenes$gene_id)

# Number of sQTLs per tissue
count(sqtls, tissue, sort = TRUE)

# sGenes found in N tissues:
sgenes |>
    summarise(n_tissues = n(), .by = gene_id) |>
    count(n_tissues) |>
    mutate(frac = n / sum(n))

overlap <- bind_rows(
    left_join(
        egenes |>
            group_by(tissue) |>
            summarise(egenes = list(gene_id)),
        sgenes |>
            group_by(tissue) |>
            summarise(sgenes = list(gene_id)),
        by = "tissue"
    ),
    left_join(
        egenes |>
            mutate(tissue = "All") |>
            group_by(tissue) |>
            summarise(egenes = list(unique(gene_id))),
        sgenes |>
            mutate(tissue = "All") |>
            group_by(tissue) |>
            summarise(sgenes = list(unique(gene_id))),
        by = "tissue"
    )
) |>
    rowwise() |>
    mutate(e_only = sum(!(egenes %in% sgenes)),
           both = sum(egenes %in% sgenes),
           s_only = sum(!(sgenes %in% egenes)),
           s_only_frac = s_only / length(sgenes)) |>
    ungroup()
overlap
overlap |>
    filter(tissue != "All") |>
    with(mean(s_only_frac))

########################################################################################
## How many expressed genes showed alternative splicing, and what fraction had sQTLs? ##
########################################################################################

# These tables are already filtered, so I'll just count the genes in them:
expr <- sgenes |>
    distinct(tissue) |>
    reframe(
        read_tsv(
            str_glue("data/expression/ensembl-gene_log2_{tissue}.bed.gz"),
            col_types = cols(gene_id = 'c', .default = '-')
        ),
        .by = tissue
    )

sj <- sgenes |>
    distinct(tissue) |>
    reframe(
        read_tsv(
            str_glue("data/splice/{tissue}.leafcutter.phenotype_groups.txt"),
            col_names = c("junction", "gene_id"),
            col_types = "cc"
        ),
        .by = tissue
    ) |>
    distinct(tissue, gene_id)

expr |>
    distinct(tissue, gene_id) |>
    left_join(
        mutate(sj, alt = TRUE), by = c("tissue", "gene_id")
    ) |>
    left_join(
        mutate(sgenes, sqtl = TRUE), by = c("tissue", "gene_id")
    ) |>
    summarise(n_expr = n(),
              n_alt = sum(!is.na(alt)),
              n_sgene = sum(!is.na(sqtl)),
              frac_sgene = n_sgene / n_alt,
              .by = tissue)
