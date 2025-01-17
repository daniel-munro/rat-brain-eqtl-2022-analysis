library(tidyverse)
library(patchwork)

d_afc <- read_tsv("analysis/orthologs/ortholog_aFC.tsv", col_types = "ccdd")
d_h2 <- read_tsv("analysis/orthologs/ortholog_h2.tsv", col_types = "ccdd")
d_sdg <- read_tsv("analysis/orthologs/ortholog_SDg.tsv", col_types = "ccdd")
d <- bind_rows(
    d_afc |>
        rename(rat = mean_abs_rat,
               human = mean_abs_human) |>
        mutate(stat = "aFC"),
    d_h2 |>
        rename(rat = mean_h2_rat,
               human = mean_h2_human) |>
        mutate(stat = "h2"),
    d_sdg |>
        rename(rat = SDg_rat,
               human = SDg_GTEx) |>
        mutate(stat = "SDg"),
)

human_genes <- read_tsv("data/gtex/gtex_genes.txt", col_types = "cc", col_names = c("id", "name")) |>
    select(name, id) |>
    deframe()

gsets <- read_tsv("data/gtex/GeneSets.txt", col_types = cols(.default = "c"), comment = "#") |>
    mutate(across(-starts_with("Autosomal"), ~ human_genes[.x])) |>
    pivot_longer(everything(), names_to = "set", values_to = "gene_id_human") |>
    filter(!is.na(gene_id_human)) |>
    mutate(set = str_replace_all(set, "zz", "")) |>
    distinct() |>
    bind_rows(tibble(set = "All (unfiltered)",
                     gene_id_human = unique(d$gene_id_human)))

gsets |>
    group_by(set) |>
    summarise(n_in_eQTL_orthologs = sum(gene_id_human %in% d$gene_id_human),
              n_total = n(),
              .groups = "drop") |>
    mutate(frac = n_in_eQTL_orthologs / n_total) |>
    arrange(frac) |>
    print.data.frame()

d_sets <- d |>
    inner_join(gsets, by = "gene_id_human")

stats <- d_sets |>
    group_by(set, stat) |>
    summarise(
        pairs = n(),
        pearson = cor(rat, human),
        pearson_p = cor.test(rat, human)$p.value,
        .groups = "drop"
    ) |>
    rowwise() |>
    mutate(
        lab = str_c(str_glue(" r = {format(pearson, digits = 2)}"),
                    str_glue("P = {format(pearson_p, digits = 2)}"),
                    str_glue("n = {pairs}"),
                    sep = "\n")
    ) |>
    ungroup()

# d_sets |>
#     ggplot(aes(x = rat, y = human)) +
#     facet_grid(rows = vars(set), cols = vars(stat), scales = "free_x") +
#     geom_point(size = 1, alpha = 0.5) +
#     coord_fixed() +
#     theme_bw()

p1 <- d_sets |>
    filter(stat == "aFC") |>
    ggplot(aes(x = rat, y = human)) +
    facet_grid(rows = vars(set), cols = vars(stat)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_text(aes(label = lab), data = filter(stats, stat == "aFC"),
              x = 3.5, y = 4.5, hjust = 0) +
    coord_fixed() +
    theme_bw()
p2 <- d_sets |>
    filter(stat == "h2") |>
    ggplot(aes(x = rat, y = human)) +
    facet_grid(rows = vars(set), cols = vars(stat)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_text(aes(label = lab), data = filter(stats, stat == "h2"),
              x = 0.6, y = 0.7, hjust = 0) +
    coord_fixed() +
    theme_bw()
p3 <- d_sets |>
    filter(stat == "SDg") |>
    ggplot(aes(x = rat, y = human)) +
    facet_grid(rows = vars(set), cols = vars(stat)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_text(aes(label = lab), data = filter(stats, stat == "SDg"),
              x = 0.3, y = 0.4, hjust = 0) +
    coord_fixed() +
    theme_bw()

p1 + p2 + p3
ggsave("analysis/orthologs/ortholog_gene_sets.png", width = 10, height = 50, dpi = 100, limitsize = FALSE)
 