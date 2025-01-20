library(argparser)
suppressPackageStartupMessages(library(tidyverse))
library(patchwork)

p <- arg_parser("Plot haplotype probabilities across chromosomes.")
p <- add_argument(p, "infile", help="name of RDS input file")
p <- add_argument(p, "outfile", help="name of PNG output file")
p <- add_argument(p, "--title", help="title to display in the plot", default=NULL)
p <- add_argument(p, "--cluster", help="whether to arrange rows with dendrogram",
                  default=TRUE)
argv <- parse_args(p)

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
    d <- cubelyr::as.tbl_cube(prob) |>
        as_tibble() |>
        mutate(
            individual = as.integer(factor(individual, levels = individuals)),
            strain = factor(strains[strain], levels = strains),
            SNP = as.integer(factor(SNP, levels = SNPs))
        )
    d
}

individuals_df <- function(d) {
    # "individual" as row, strain_SNP as columns.
    df <- d |>
        pivot_wider(names_from = c(strain, SNP), values_from = prob) |>
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
    )
}

apr <- readRDS(argv$infile)
d <- probs(apr)

if (argv$cluster) {
    clust <- hclust(dist(individuals_df(d)))
    ddata <- ggdendro::dendro_data(as.dendrogram(clust))
    p1 <- ggplot(ggdendro::segment(ddata)) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        coord_flip() +
        scale_y_reverse() +
        scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
        ggdendro::theme_dendro()

    p2 <- d |>
        mutate(individual = factor(individual, levels = rev(ggdendro::label(ddata)$label))) |>
        plot_probs() +
        ylab(NULL) +
        ggtitle(argv$title)
    
    p1 + p2 + plot_layout(widths = c(1, 10))
} else {
    plot_probs(d) +
        ylab(NULL) +
        ggtitle(argv$title)
}

height <- if (dim(apr)[1] > 60) 24 else 12
ggsave(argv$outfile, width = 8, height = height)
