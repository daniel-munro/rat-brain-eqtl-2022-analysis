library(VariantAnnotation)
library(tidyverse)

r2 <- function(m1, m2) {
    # for each row i, computes correlation between m1[i, ] and m2[i, ]
    a1 <- m1 - rowMeans(m1)
    a2 <- m2 - rowMeans(m2)
    (rowSums(a1 * a2) / sqrt(rowSums(a1 ^ 2) * rowSums(a2 ^ 2))) ^ 2
}

esnps <- read_tsv("data/eqtls/top_assoc.txt", col_types = "ccciiciiccdddddddd") |>
    filter(qval < 0.05) |>
    select(tissue, gene_id, variant_id)

filename <- "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
gt <- VariantAnnotation::readGT(filename)
gt <- gt[unique(esnps$variant_id), ]
geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)[x])
rownames(geno) <- rownames(gt)
rm(gt)


ld_pairs <- esnps |>
    group_by(gene_id) |>
    summarise(
        crossing(tissue.x = unique(tissue),
                 tissue.y = unique(tissue)) |>
            filter(tissue.x < tissue.y) |>
            left_join(tibble(tissue.x = tissue,
                             variant_id.x = variant_id),
                      by = "tissue.x") |>
            left_join(tibble(tissue.y = tissue,
                             variant_id.y = variant_id),
                      by = "tissue.y")
    )

ld <- crossing(gene_id = unique(esnps$gene_id),
                     tissue.x = unique(esnps$tissue),
                     tissue.y = unique(esnps$tissue)) |>
    filter(tissue.x < tissue.y) |>
    left_join(rename(esnps, tissue.x = tissue),
              by = c("gene_id", "tissue.x")) |>
    left_join(rename(esnps, tissue.y = tissue),
              by = c("gene_id", "tissue.y")) |>
    filter(!is.na(variant_id.x),
           !is.na(variant_id.y)) |>
    mutate(r2 = r2(geno[variant_id.x, , drop = FALSE],
                   geno[variant_id.y, , drop = FALSE]))

# write_tsv(ld, "specificity/tissue_eSNP_LD.txt")

ld |>
    ggplot(aes(x = r2)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = median(ld$r2)) +
    geom_vline(xintercept = mean(ld$r2), color = "red") +
    ggtitle("LD between top eSNPs for same gene in different tissues",
            subtitle = str_glue("Median: {median(ld$r2) |> round(2)}, Mean: {mean(ld$r2) |> round(2)}, {round(100 * mean(ld$r2 == 1))}% are 1, {round(100 * mean(ld$r2 > 0.6))}% are above 0.6"))

ggsave("specificity/tissue_eSNP_LD.png", width = 6, height = 4)
