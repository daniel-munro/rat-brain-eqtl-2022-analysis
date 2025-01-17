library(tidyverse)
library(patchwork)

probs <- function(prob) {
    names(dimnames(prob)) <- c("individual", "Strain", "SNP")
    cubelyr::as.tbl_cube(prob) |>
        as_tibble() |>
        mutate(Strain = strains[Strain])
}

# plot_mean_probs <- function(chr) {
#     mean_probs |>
#         filter(chrom == chr) |>
#         ggplot(aes(x = pos, y = prob, fill = Strain)) +
#         geom_area() +
#         scale_fill_brewer(type = "qual", palette = 6) +
#         theme_minimal() +
#         scale_x_continuous(expand = c(0, 0)) +
#         theme(
#             axis.text.y = element_blank(),
#             panel.grid.major.y = element_blank(),
#             panel.grid.minor.y = element_blank(),
#         ) +
#         xlab(NULL) +
#         ylab(NULL)
# }
# 
# plot_uncert <- function(chr) {
#     uncert |>
#         filter(chrom == chr) |>
#         ggplot(aes(x = pos, y = entropy)) +
#         geom_smooth(
#             stat = "summary",
#             fun.data = function(yval) tibble(ymin = quantile(yval, 0.05),
#                                              y = median(yval),
#                                              ymax = quantile(yval, 0.95)),
#             size = 0.5
#         ) +
#         theme_minimal() +
#         scale_x_continuous(expand = c(0, 0)) +
#         xlab(NULL) +
#         ylab(NULL)
# }

strains <- c(
    A = "ACI/N", B = "BN/SsN", C = "BUF/N", D = "F344/N",
    E = "M520/N", F = "MR/N", G = "WN/N", H = "WKY/N"
)

d <- tibble(chrom = 1:20) |>
    group_by(chrom) |>
    summarise(
        readRDS(str_glue("data/haplotypes/haplotype_probs_chr{chrom}.rds")) |>
            probs(),
        .groups = "drop"
    )

mean_probs <- d |>
    group_by(SNP, Strain) |>
    summarise(prob = mean(prob), .groups = "drop") |>
    separate(SNP, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    # filter(pos < 3e7) |>
    mutate(pos = pos / 1e6,
           chrom = factor(chrom, levels = str_c("chr", 1:20)))

uncert <- d |>
    group_by(SNP, individual) |>
    summarise(entropy = entropy::entropy.empirical(prob, unit = "log2"),
              .groups = "drop") |>
    separate(SNP, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    # filter(pos < 3e7) |>
    mutate(pos = pos / 1e6,
           chrom = factor(chrom, levels = str_c("chr", 1:20)))

mean_probs |>
    # filter(chrom %in% str_c("chr", 1:6)) |>
    ggplot(aes(x = pos, y = prob, fill = Strain)) +
    facet_wrap(~ chrom, ncol = 2, dir = "v", strip.position = "left") +
    geom_area() +
    scale_fill_brewer(type = "qual", palette = 6) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    # xlab("Position (Mb)") +
    # ylab("Mean probability") +
    xlab(NULL) +
    ylab(NULL) +
    theme(
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing.y = unit(2.5, "lines"),
        legend.position = "top",
    )

ggsave("haplotypes/paper/hap_combined_1.png", width = 11, height = 9)

uncert |>
    # filter(chrom %in% str_c("chr", 1:6)) |>
    ggplot(aes(x = pos, y = entropy)) +
    facet_wrap(~ chrom, ncol = 2, dir = "v", strip.position = "left") +#, labeller = labeller(chrom = "")) +
    geom_smooth(
        stat = "summary",
        fun.data = function(yval) tibble(ymin = quantile(yval, 0.05),
                                         y = median(yval),
                                         ymax = quantile(yval, 0.95)),
        size = 0.5
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    xlab("Position (Mb)") +
    # ylab("Entropy (base 2) per animal") +
    # xlab(NULL) +
    ylab(NULL) +
    theme_minimal() +
    theme(
        panel.spacing.y = unit(2.5, "lines"),
        # strip.text = element_blank(),
    )

ggsave("haplotypes/paper/hap_combined_2.png", width = 11, height = 8.25)

# p1 <- map(str_c("chr", 1:3), plot_mean_probs)
# p2 <- map(str_c("chr", 1:3), plot_uncert)

# p1[[1]] / p2[[1]] / p1[[2]] / p2[[2]] / p1[[3]] / p2[[3]]
# p1[[1]] / p2[[1]]
