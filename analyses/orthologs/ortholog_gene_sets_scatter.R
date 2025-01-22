library(tidyverse)
library(patchwork)

d_afc <- read_tsv("data/gtex/orthologs/ortholog_aFC.tsv", col_types = "ccdd")
d_h2 <- read_tsv("data/gtex/orthologs/ortholog_h2.tsv", col_types = "ccdd")
d_sdg <- read_tsv("data/gtex/orthologs/ortholog_SDg.tsv", col_types = "ccdd")
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

gsets <- read_tsv("data/gtex/orthologs/gene_sets.tsv", col_types = "cc") |>
    bind_rows(tibble(set = "All (unfiltered)",
                     gene_id_human = unique(d$gene_id_human)))

gsets |>
    summarise(n_in_eQTL_orthologs = sum(gene_id_human %in% d$gene_id_human),
              n_total = n(),
              .by = set) |>
    mutate(frac = n_in_eQTL_orthologs / n_total) |>
    arrange(frac) |>
    print.data.frame()

d_sets <- d |>
    inner_join(gsets, by = "gene_id_human", relationship = "many-to-many")

stats <- d_sets |>
    summarise(
        pairs = n(),
        pearson = cor(rat, human),
        pearson_p = cor.test(rat, human)$p.value,
        .by = c(set, stat)
    ) |>
    rowwise() |>
    mutate(
        lab = str_c(str_glue(" r = {format(pearson, digits = 2)}"),
                    str_glue("P = {format(pearson_p, digits = 2)}"),
                    str_glue("n = {pairs}"),
                    sep = "\n")
    ) |>
    ungroup()

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
ggsave("analyses/orthologs/ortholog_gene_sets_scatter.png", width = 10, height = 50, dpi = 100, limitsize = FALSE)
