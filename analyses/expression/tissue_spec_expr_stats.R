library(tidyverse)

samples <- read_tsv("data/samples.txt", col_types = "cc") |>
    rename(tissue = brain_region)

expr <- read.table("data/expression/ensembl-gene_raw-counts.txt",
                   check.names = FALSE) |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "library", values_to = "expr") |>
    filter(library %in% samples$library) |>
    left_join(samples, by = "library") |>
    summarise(median_expr = median(expr),
              frac_nonzero = mean(expr > 0),
              frac_expr = mean(expr >= 10),
              .by = c(tissue, gene_id))

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    distinct(tissue, gene_id)

expr |>
    filter(gene_id %in% eqtls$gene_id) |>
    left_join(mutate(eqtls, eqtl = TRUE), by = c("tissue", "gene_id")) |>
    replace_na(list(eqtl = FALSE)) |>
    summarise(mean_frac_nonzero = mean(frac_nonzero),
              mean_frac_expr = mean(frac_expr),
              .by = eqtl)
