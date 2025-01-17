suppressPackageStartupMessages(library(tidyverse))

tissues <- c("Acbc", "IL", "LHB", "PL", "VoLo")

expr <- tibble(tissue = tissues) %>%
    group_by(tissue) %>%
    summarise(
        read_tsv(str_glue("data/expression/ensembl-gene_log2_{tissue}.bed.gz"),
                 col_types = cols(`#chr` = "-", start = "-", end = "-",
                                  gene_id = "c", .default = "d")) %>%
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "expr"),
        .groups = "drop"
    ) %>%
    # Convert from log2(count + 1) to count.
    mutate(expr = round(2 ^ expr - 1)) %>%
    group_by(tissue, gene_id) %>%
    summarise(median_expr = median(expr),
              # expressed = median_expr > 0,
              frac_nonzero = mean(expr > 0),
              # expressed = frac_nonzero > 0.75,
              frac_expr = mean(expr >= 10),
              .groups = "drop") %>%
    # Add back filtered genes with maximum possible expression given the filter,
    # making sure the particular values don't affect stats:
    mutate(present = TRUE) %>%
    complete(tissue, gene_id, fill = list(present = FALSE, frac_nonzero = 0.25, frac_expr = 0.25))

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") %>%
    distinct(tissue, gene_id)

# genes <- eqtls %>%
#     # left_join(expr, by = c("tissue", "gene_id")) %>%
#     group_by(gene_id) %>%
#     summarise(n_tissues_eqtl = n_distinct(tissue),
#               # n_expressed = sum(expressed),
#               .groups = "drop")

# counts <- genes %>%
#     count(n_tissues, n_expressed)

expr %>%
    filter(gene_id %in% eqtls$gene_id) %>%
    # left_join(genes, by = "gene_id") %>%
    left_join(mutate(eqtls, eqtl = TRUE), by = c("tissue", "gene_id")) %>%
    replace_na(list(eqtl = FALSE)) %>%
    group_by(eqtl) %>%
    summarise(mean_frac_nonzero = mean(frac_nonzero),
              mean_frac_expr = mean(frac_expr),
              .groups = "drop")
    
    # filter(n_tissues_eqtl == 1) %>%
    # group_by(gene_id) %>%
    # summarise(n_tissues_med_expr = sum(frac_expr > 0.5),
    #           .groups = "drop")
# tmp %>%
#     count(n_tissues_med_expr)
# tmp %>%
#     summarise(frac_all_expr = mean(n_tissues_med_expr == 5),
#               .groups = "drop")
