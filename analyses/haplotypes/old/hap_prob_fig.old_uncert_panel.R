# Based on ../haplotype_skew.R

library(tidyverse)

probs <- function(prob) {
    names(dimnames(prob)) <- c("individual", "Strain", "SNP")
    cubelyr::as.tbl_cube(prob) |>
        as_tibble() |>
        mutate(Strain = strains[Strain])
}

strains <- c(
    A = "ACI/N", B = "BN/SsN", C = "BUF/N", D = "F344/N",
    E = "M520/N", F = "MR/N", G = "WN/N", H = "WKY/N"
)
# chr_len <- c(282763074, 266435125, 177699992)
# n_SNPs <- round(2000 * chr_len / chr_len[1])

d <- tibble(chrom = 1:3) |>
    group_by(chrom) |>
    summarise(
        readRDS(str_glue("data/haplotypes/haplotype_probs_chr{chrom}.rds")) |>
            probs(),
        .groups = "drop"
    ) |>
    group_by(SNP, Strain) |>
    summarise(prob = mean(prob), .groups = "drop") |>
    separate(SNP, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(pos = pos / 1e6)

ggplot(d, aes(x = pos, y = prob, fill = Strain)) +
    facet_wrap(~ chrom, ncol = 1, strip.position = "left") +
    geom_area() +
    scale_fill_brewer(type = "qual", palette = 6) +
    theme_minimal() +
    scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
    ) +
    xlab("Position (Mb)") +
    ylab("Mean probability")

ggsave("haplotypes/paper/hap_prob.png", width = 6, height = 2.5)

############################################
## Uncertainty (individual-level entropy) ##
############################################

d2 <- tibble(chrom = 1:3) |>
    group_by(chrom) |>
    summarise(
        readRDS(str_glue("data/haplotypes/haplotype_probs_chr{chrom}.rds")) |>
            probs(),
        .groups = "drop"
    ) |>
    group_by(SNP, individual) |>
    summarise(entropy = entropy::entropy.empirical(prob, unit = "log2"),
              .groups = "drop") |>
    # group_by(SNP) |>
    # summarise(entropy_05 = quantile(entropy, 0.05),
    #           entropy_50 = median(entropy),
    #           entropy_95 = quantile(entropy, 0.95)) |>
    separate(SNP, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(pos = pos / 1e6)

ggplot(d2, aes(x = pos, y = entropy)) +
    facet_wrap(~ chrom, ncol = 1, strip.position = "right") +
    # facet_wrap(~ chrom, ncol = 1) +
    geom_smooth(
        stat = "summary",
        fun.data = function(yval) tibble(ymin = quantile(yval, 0.05),
                                         y = median(yval),
                                         ymax = quantile(yval, 0.95)),
        size = 0.5
    ) +
    theme_minimal() +
    scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0)) +
    xlab("Position (Mb)") +
    ylab("Entropy (base 2) per animal")

ggsave("haplotypes/paper/hap_uncert.png", width = 5, height = 2.5)
