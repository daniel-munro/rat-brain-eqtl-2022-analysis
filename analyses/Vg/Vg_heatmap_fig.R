library(tidyverse)

heatmap_order <- function(fac1, fac2, value) {
    df <- tibble(fac1, fac2, value) |>
        pivot_wider(fac1, names_from = fac2, values_from = value)
    mat <- df |> select(-fac1) |> as.matrix()
    rownames(mat) <- df$fac1
    clust <- hclust(dist(mat))
    factor(fac1, levels = clust$labels[clust$order])
}

gtex <- read_tsv("analysis/Vg/gtex_phaser_vg.tsv.gz",
                 col_types = cols(GENE_ID = "c", .default = "d")) |>
    rename(gene_id = GENE_ID) |>
    pivot_longer(-gene_id, names_to = "tissue", values_to = "Vg") |>
    filter(str_sub(tissue, 1, 3) == "BRN",
           !is.na(Vg))

rat <- tibble(tissue = c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("analysis/Vg/anevah/Vg.{tissue}.tsv.gz"), col_types = "cd"),
        .groups = "drop"
    ) |>
    rename(gene_id = GENE_ID)

gene_map <- read_tsv("data/gtex/orthologs.txt", col_types = "cc") |>
    mutate(gene_id_human = str_replace(gene_id_human, "\\..+$", "")) |>
    filter(gene_id_rat %in% rat$gene_id,
           gene_id_human %in% gtex$gene_id)

vg <- bind_rows(gtex, rat) |>
    filter(gene_id %in% c(gene_map$gene_id_rat, gene_map$gene_id_human))

tissues <- jsonlite::fromJSON("data/gtex/tissueInfo.json")[[1]] |>
    select(tissueSiteDetailAbbr,
           tissueSiteDetail) |>
    mutate(tissueSiteDetail = tissueSiteDetail |>
               str_replace("Brain - ", "") |>
               str_c(" - Human")) |>
    deframe() |>
    c(IL = "Infralimbic cortex - Rat",
      LHb = "Lateral habenula - Rat",
      NAcc = "Nucleus accumbens core - Rat",
      OFC = "Orbitofrontal cortex - Rat",
      PL = "Prelimbic cortex - Rat")

cor_intra <- bind_rows(
    crossing(tissue.x = unique(rat$tissue),
             tissue.y = unique(rat$tissue)),
    crossing(tissue.x = unique(gtex$tissue),
             tissue.y = unique(gtex$tissue))
) |>
    group_by(tissue.x, tissue.y) |>
    summarise(
        inner_join(
            filter(vg, tissue == tissue.x),
            filter(vg, tissue == tissue.y),
            by = "gene_id"
        ) |>
            summarise(r = cor(Vg.x, Vg.y),
                      n = n()),
        .groups = "drop"
    )

cor_inter <- crossing(tissue.x = unique(rat$tissue),
                      tissue.y = unique(gtex$tissue)) |>
    group_by(tissue.x, tissue.y) |>
    summarise(
        gene_map |>
            inner_join(
                filter(vg, tissue == tissue.x),
                by = c("gene_id_rat" = "gene_id")
            ) |>
            inner_join(
                filter(vg, tissue == tissue.y),
                by = c("gene_id_human" = "gene_id")
            ) |>
            summarise(r = cor(Vg.x, Vg.y),
                      n = n()),
        .groups = "drop"
    )

cor_inter <- bind_rows(
    cor_inter,
    cor_inter |>
        mutate(tmp = tissue.x) |>
        mutate(tissue.x = tissue.y,
               tissue.y = tmp) |>
        select(-tmp)
)

bind_rows(cor_intra, cor_inter) |>
    mutate(tissue.x = tissues[tissue.x],
           tissue.y = tissues[tissue.y],
           tissue.x = heatmap_order(tissue.x, tissue.y, r),
           tissue.y = heatmap_order(tissue.y, tissue.x, r)) |>
    filter(as.integer(tissue.x) <= as.integer(tissue.y)) |>
    mutate(tissue.y = fct_rev(tissue.y)) |>
    ggplot(aes(x = tissue.x, y = tissue.y, fill = r)) +
    geom_tile() +
    coord_fixed() +
    expand_limits(fill = 0) +
    scale_fill_viridis_c() +
    xlab(NULL) +
    ylab(NULL) +
    labs(fill = expression(Cor(V^G))) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        panel.grid = element_blank(),
        legend.position = c(0.85, 0.7),
        legend.key.size = unit(15, "pt"),
    )

ggsave("analysis/Vg/Vg_heatmap_fig.png", width = 5, height = 5, bg = "white")

###########
## Stats ##
###########

sd(rat$Vg)
sd(gtex$Vg)
median(rat$Vg)
median(gtex$Vg)
summary(sqrt(rat$Vg))
summary(sqrt(gtex$Vg))
