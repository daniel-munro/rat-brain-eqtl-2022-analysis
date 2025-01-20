library(VariantAnnotation)
library(tidyverse)

################################
## LD plot for example region ## (updated version below)
################################

load_geno <- function(chrom, start, end) {
    filename <- "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
    rng <- GRanges(chrom, IRanges(start, end))
    gt <- readGT(filename, param = rng)
    geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)[x])
    rownames(geno) <- rownames(gt)
    t(geno)
}

## LD across whole chromosome

geno <- load_geno("12", 1, 536870912)     # For every nth chr12 snp
all_snps <- colnames(geno)
geno <- geno[, apply(geno, 2, var) != 0]
geno <- geno[, seq(1, ncol(geno), by = ncol(geno) %/% 1000)]
ld <- cor(geno) ^ 2

d_ld <- ld |>
    as_tibble(rownames = "var.x") |>
    pivot_longer(-var.x, names_to = "var.y", values_to = "LD") |>
    separate(var.x, c("chr.x", "pos.x"), convert = TRUE) |>
    separate(var.y, c("chr.y", "pos.y"), convert = TRUE)

d_ld |>
    filter(pos.x <= pos.y) |>
    mutate(LD = pmax(0.2, LD)) |>
    ggplot(aes(x = as.factor(pos.x), y = as.factor(pos.y), fill = LD)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black", limits = c(0.2, 1)) +
    coord_fixed() +
    xlab(NULL) +
    ylab(NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank())

ggsave("analyses/genotypes/LD_top.png", width = 4, height = 4, dpi = 300)

## Links from whole-chromosome heatmap to chromosome

tibble(pos = as.integer(str_split(all_snps, ":", simplify = TRUE)[, 2])) |>
    mutate(rank = rank(pos),
           pos = pos / 1e6) |>
    ggplot(aes(x = pos, xend = rank * (max(pos) / max(rank)), y = 0, yend = 1)) +
    geom_segment(size = 0.01, color = "#aaaacc") +
    ylab(NULL) +
    xlab(NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank())

ggsave("analyses/genotypes/LD_middle.png", width = 4, height = 1, dpi = 300, bg = "white")

## LD for first 1000 SNPs on chromosome or first 2 Mb

geno2 <- load_geno("12", 1, 2e6)
ld2 <- cor(geno2) ^ 2

d_ld2 <- ld2 |>
    as_tibble(rownames = "var.x") |>
    pivot_longer(-var.x, names_to = "var.y", values_to = "LD") |>
    separate(var.x, c("chr.x", "pos.x"), convert = TRUE) |>
    separate(var.y, c("chr.y", "pos.y"), convert = TRUE)

d_ld2 |>
    filter(pos.x >= pos.y) |>
    ggplot(aes(x = as.factor(pos.x), y = as.factor(pos.y), fill = LD)) +
    geom_tile() +
    scale_fill_viridis_c(direction = -1, option = "A", limits = c(-0.001, 1)) +
    coord_fixed() +
    xlab(NULL) +
    ylab(NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank())

ggsave("analyses/genotypes/LD_bottom.png", width = 4, height = 4, dpi = 300)

## Links from segment heatmap to chromosome

tibble(pos = unique(d_ld2$pos.x)) |>
    mutate(rank = rank(pos),
           pos = pos / 1e6) |>
    ggplot(aes(x = pos, xend = rank * (max(d_ld$pos.x / 1e6) / max(rank)), y = 1, yend = 0)) +
    geom_segment(size = 0.05, color = "#aaaacc") +
    ylab(NULL) +
    xlab(NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank())

ggsave("analyses/genotypes/LD_middle2.png", width = 4, height = 0.7, dpi = 300, bg = "white")

##########################################
## Simplified LD plot: first cis-window ##
##########################################

tss <- read_table("data/genes.txt", col_types = "c-c---illlll") |>
    filter(in_expr_IL | in_expr_LHb | in_expr_NAcc | in_expr_OFC | in_expr_PL,
           chrom == "1") |>
    filter(tss == min(tss)) |>
    pull(tss)

geno <- load_geno("1", tss - 1e6, tss + 1e6)
ld <- cor(geno) ^ 2

d_ld <- ld |>
    as_tibble(rownames = "var.x") |>
    pivot_longer(-var.x, names_to = "var.y", values_to = "r2") |>
    mutate(pos.x = var.x |> str_replace("chr1:", "") |> as.integer(),
           pos.y = var.y |> str_replace("chr1:", "") |> as.integer())

d_ld |>
    filter(pos.x >= pos.y) |>
    ggplot(aes(x = as.factor(pos.x), y = as.factor(pos.y), fill = r2)) +
    geom_tile() +
    scale_fill_viridis_c(direction = -1, option = "A", limits = c(-0.001, 1)) +
    coord_fixed() +
    xlab(NULL) +
    ylab(NULL) +
    labs(fill = expression("LD "*(r^2))) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank())

ggsave("analyses/genotypes/LD_ciswindow.png", width = 4, height = 4, dpi = 300)
