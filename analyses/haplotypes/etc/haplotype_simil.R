suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
library(reticulate)

args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]

strains <- c(
    A = "ACI",
    B = "BN",
    C = "BUF",
    D = "F344",
    E = "M520",
    F = "MR",
    G = "WN",
    H = "WKY"
)

# pair_index_to_probs <- function(pair_index) {
#     apr <- array(0, dim = c(dim(pair_index)[1], 8, dim(pair_index)[2]))
#     index_pairs <- c()
#     for (i in 1:7) {
#         for (j in ((i + 1):8)) {
#             index_pairs <- c(index_pairs
#         }
#     }
# }

probs <- function(prob, n_SNPs = NULL, n_individuals = NULL) {
    names(dimnames(prob)) <- c("individual", "strain", "SNP")
    if (!is.null(n_individuals)) {
        individuals <- sort(sample(dim(apr)[1], n_individuals, replace = FALSE))
        prob <- prob[individuals, , ]
    }
    if (!is.null(n_SNPs)) {
        SNP_reps <- round(seq(from = 1, to = dim(prob)[3], length.out = n_SNPs))
        prob <- prob[, , SNP_reps]
    }
    individuals <- dimnames(prob)[["individual"]]
    SNPs <- dimnames(prob)[["SNP"]]
    d <- as.tbl_cube(prob) %>%
        as_tibble() %>%
        mutate(
            individual = as.integer(factor(individual, levels = individuals)),
            strain = factor(strains[strain], levels = strains),
            SNP = as.integer(factor(SNP, levels = SNPs))
        )
    d
}

order_individuals <- function(d) {
    df <- individuals_df(d)
    pca <- prcomp(df)
    rownames(df)[order(pca$x[, "PC1"])]
}

individuals_df <- function(d) {
    # "individual" as row, strain_SNP as columns.
    df <- d %>%
        pivot_wider(names_from = c(strain, SNP), values_from = prob) %>%
        as.data.frame()
    rownames(df) <- df$individual
    df$individual <- NULL
    df
}

plot_probs <- function(d) {
    ggplot(d, aes(x = SNP, y = prob, fill = strain)) +
        facet_grid(rows = "individual") +
        geom_col(width = 1) +
        scale_fill_brewer(type = "qual", palette = 6) +
        theme_minimal() +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(
            strip.text = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank(),
            axis.text.x = element_blank(),
            # panel.spacing = unit(0.002, "npc")
        )
}

np <- import("numpy")
apr <- np$load(str_c("haplotypes/haplotype_probs/haplotype_sim_chr", chrom, ".npy"))
# apr <- pair_index_to_apr(apr)
dimnames(apr) <- list(1:dim(apr)[1], names(strains), 1:dim(apr)[3])
d <- probs(apr)

clust <- hclust(dist(individuals_df(d)))
ddata <- ggdendro::dendro_data(as.dendrogram(clust))
p1 <- ggplot(ggdendro::segment(ddata)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_flip() +
    scale_y_reverse() +
    scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
    ggdendro::theme_dendro()
# geom_text(data = ggdendro::label(ddata), 
#           aes(x = x, y = y, label = label), vjust = -0.5, size = 3)

p2 <- d %>%
    mutate(individual = factor(individual, levels = rev(ggdendro::label(ddata)$label))) %>%
    plot_probs() +
    ylab(NULL) +
    ggtitle(str_c("Chromosome ", chrom))

p1 + p2 + plot_layout(widths = c(1, 10))

ggsave(str_c("haplotypes/chr_plots/sim_chr", chrom, ".png"), width = 8, height = 24)
