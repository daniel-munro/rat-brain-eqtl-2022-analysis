library(tidyverse)

rat <- read_tsv("data/expr_var_explained/expr_R2.tsv", col_types = "ccdd") |>
    select(-corr)

human <- read_tsv(list.files("data/expr_var_explained/R2_brain", full.names = TRUE),
                  col_types = "cd", id = "path") |>
    mutate(tissue = str_match(path, ".*Brain_([^/]+)_R2.tsv")[, 2]) |>
    select(tissue,
           gene_id = Name,
           R2 = R2_all_variants) |>
    filter(!is.na(R2))

# Made in rat_human_aFC_stats.R:
orthologs <- read_tsv("data/afc/gene_map.txt", col_types = "cc") |>
    mutate(gene_id_human = str_replace(gene_id_human, "\\..+$", ""))

orthologs_filt <- orthologs |>
    inner_join(
        count(rat, gene_id, name = "n_rat_tissues"),
        by = c("gene_id_rat" = "gene_id")
    ) |>
    inner_join(
        count(human, gene_id, name = "n_human_tissues"),
        by = c("gene_id_human" = "gene_id")
    ) |>
    filter(n_rat_tissues == n_distinct(rat$tissue),
           n_human_tissues == n_distinct(human$tissue))

r2 <- bind_rows(
    rat |>
        filter(gene_id %in% orthologs_filt$gene_id_rat) |>
        mutate(Organism = "Rat", .before = 1),
    human |>
        filter(gene_id %in% orthologs_filt$gene_id_human) |>
        mutate(Organism = "Human", .before = 1)
)

tissues <- c("Infralimbic cortex" = "IL",
             "Lateral habenula" = "LHb",
             "Nucleus accumbens core" = "NAcc",
             "Orbitofrontal cortex" = "OFC",
             "Prelimbic cortex" = "PL")
r2 |>
    mutate(tissue = tissue |>
               str_replace_all("_", " ") |>
               fct_recode(!!!tissues) |>
               fct_reorder(R2, median)) |>
    ggplot(aes(y = tissue, x = R2, fill = Organism)) +
    geom_boxplot(outlier.size = 0.1) +
    xlab(expression(R^2*" for ortholog-filtered genes")) +
    ylab(NULL) +
    theme_minimal() +
    ggtitle("Only the 66 ortholog pairs with R2 for all tissues")

ggsave("analyses/expr_var_explained/expr_R2.png", height = 3.5, width = 8, bg = "white")

r2 |>
    group_by(Organism, tissue) |>
    skimr::skim(R2)

##############################################################
## Correlation in average rat vs human R2 per ortholog pair ##
##############################################################

rat_mean <- rat |>
    group_by(gene_id) |>
    summarise(mean_R2_rat = mean(R2))

human_mean <- human |>
    group_by(gene_id) |>
    summarise(mean_R2_human = mean(R2))

R2_corr <- orthologs |>
    inner_join(rat_mean, by = c("gene_id_rat" = "gene_id")) |>
    inner_join(human_mean, by = c("gene_id_human" = "gene_id"))

corrs <- R2_corr |>
    summarise(
        R = cor(mean_R2_rat, mean_R2_human),
        rho = cor(mean_R2_rat, mean_R2_human, method = "s"),
        n = n(),
        .groups = "drop"
    ) |>
    mutate(
        stats = str_c("R = ", format(R, digits = 2, nsmall = 2),
                      "\nrho = ", format(rho, digits = 2, nsmall = 2),
                      "\nN = ", n)
    )

R2_corr |>
    ggplot(aes(x = mean_R2_rat, y = mean_R2_human)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_text(aes(x = 0.73, y = 0.87, label = stats), data = corrs, hjust = "left") +
    coord_fixed() +
    expand_limits(x = 1, y = 1) +
    xlab(expression("Mean "*R^2*" in Rat")) +
    ylab(expression("Mean "*R^2*" in Human")) +
    theme_minimal()

ggsave("analyses/expr_var_explained/expr_R2_corr.png", width = 4, height = 4, bg = "white")

with(R2_corr, cor.test(mean_R2_rat, mean_R2_human))
with(R2_corr, cor(mean_R2_rat, mean_R2_human, method = "spearman"))
