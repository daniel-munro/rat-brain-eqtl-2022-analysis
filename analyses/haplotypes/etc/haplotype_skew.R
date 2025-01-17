suppressPackageStartupMessages(library(tidyverse))

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
    # individuals <- dimnames(prob)[["individual"]]
    SNPs <- dimnames(prob)[["SNP"]]
    d <- cubelyr::as.tbl_cube(prob) %>%
        as_tibble() %>%
        mutate(
            # individual = as.integer(factor(individual, levels = individuals)),
            strain = factor(strains[strain], levels = strains)
            # SNP = factor(SNP, levels = SNPs)
        )
    d
}

d <- tibble(chrom = 1:20) %>%
    group_by(chrom) %>%
    do(
        readRDS(str_c("../data/analysis/haplotype_probs_chr", .$chrom, ".rds")) %>%
            probs() %>%
            separate(SNP, c("chr", "pos"), sep = ":", convert = TRUE) %>%
            group_by(pos, strain) %>%
            summarise(prob = mean(prob), .groups = "drop")
    ) %>%
    ungroup() %>%
    mutate(pos = as.factor(pos))

ggplot(d, aes(x = pos, y = prob, fill = strain)) +
    facet_wrap(~ chrom, ncol = 1, scales = "free_x", strip.position = "right") +
    geom_col(width = 1) +
    scale_fill_brewer(type = "qual", palette = 6) +
    theme_minimal() +
    # scale_x_discrete(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank()
    ) +
    ggtitle("Mean probabilities per SNP for all chromosomes")

ggsave("haplotypes/all_chr_probs.png", width = 8, height = 8)

geno <- read_csv("haplotypes/geno.csv", col_types = cols(.default = "c"), na = "-") %>%
    # mutate(SNP = 1:n()) %>%
    # filter(id %in% d$SNP) %>%
    separate(id, c("chrom", "pos"), sep = ":") %>%
    mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer()) %>%
    # Only include SNPs in 2000-per-chromosome sampling:
    inner_join(distinct(d, chrom, pos), by = c("chrom", "pos")) %>%
    pivot_longer(-c(chrom, pos), names_to = "individual", values_to = "genotype") %>%
    mutate(pos = factor(pos, levels = levels(d$pos)),
           # genotype = fct_relevel(genotype, "A", "H", "B"))
           genotype = c(A = 0L, H = 1L, B = 2L)[genotype])

fgeno <- read_csv("data/qtl2/founder_geno.csv", col_types = cols(.default = "c"), na = "-") %>%
    # mutate(SNP = 1:n()) %>%
    # filter(id %in% d$SNP) %>%
    separate(id, c("chrom", "pos"), sep = ":") %>%
    mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer()) %>%
    # Only include SNPs in 2000-per-chromosome sampling:
    inner_join(distinct(d, chrom, pos), by = c("chrom", "pos")) %>%
    pivot_longer(-c(chrom, pos), names_to = "strain", values_to = "genotype") %>%
    mutate(pos = factor(pos, levels = levels(d$pos)),
           # genotype = as.factor(genotype))
           # genotype = factor(genotype, levels = c("A", "H", "B")))
           genotype = c(A = 0L, H = 1L, B = 2L)[genotype])

geno %>%
    # filter(chrom <= 2) %>%
    count(chrom, pos, genotype) %>%
    complete(nesting(chrom, pos), genotype, fill = list(n = 0)) %>%
    group_by(chrom, genotype) %>%
    mutate(n = slider::slide_dbl(n, mean, .before = 10, .after = 10)) %>%
    ungroup() %>%
    group_by(chrom, pos) %>%
    mutate(fraction = n / sum(n)) %>%
    ggplot(aes(x = pos, y = fraction, fill = genotype)) +
    facet_wrap(~ chrom, ncol = 1, scales = "free_x", strip.position = "right") +
    geom_col(width = 1) +
    theme_minimal() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        # strip.text = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
    ) +
    ggtitle("Alt allele count, rolling average (21 SNP window)")

ggsave("haplotypes/all_chr_geno_rollfrac.png", width = 8, height = 8)

# Same but without rolling average, since in a region dominated by one
# haplotype, SNPs could be all 0 or all 1"
geno %>%
    # filter(chrom <= 2) %>%
    count(chrom, pos, genotype) %>%
    complete(nesting(chrom, pos), genotype, fill = list(n = 0)) %>%
    # group_by(chrom, genotype) %>%
    # mutate(n = slider::slide_dbl(n, mean, .before = 10, .after = 10)) %>%
    # ungroup() %>%
    group_by(chrom, pos) %>%
    mutate(fraction = n / sum(n)) %>%
    ggplot(aes(x = pos, y = fraction, fill = genotype)) +
    facet_wrap(~ chrom, ncol = 1, scales = "free_x", strip.position = "right") +
    geom_col(width = 1) +
    theme_minimal() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        # strip.text = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
    ) +
    ggtitle("Alt allele counts")

ggsave("haplotypes/all_chr_geno_fraction.png", width = 8, height = 8)

fgeno %>%
    # filter(chrom <= 2) %>%
    filter(genotype == 2) %>%
    ggplot(aes(x = pos, y = fct_rev(strain), fill = strain)) +
    facet_wrap(~ chrom, ncol = 1, scales = "free_x", strip.position = "right") +
    geom_tile() +
    scale_fill_brewer(type = "qual", palette = 6) +
    # scale_fill_discrete(drop = FALSE) +
    theme_minimal() +
    # scale_y_discrete(limits = levels(fgeno$strain)) +
    theme(
        # strip.text = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
    ) +
    ylab("Strain")

ggsave("haplotypes/all_chr_fgeno.png", width = 8, height = 8)

# At each position, what fraction of sample genotypes match each of 36 founder pairs?

strain_pairs <- crossing(strain1 = unique(fgeno$strain),
                         strain2 = unique(fgeno$strain)) %>%
    filter(strain1 <= strain2)

fpairs <- strain_pairs %>%
    left_join(fgeno, by = c("strain1" = "strain")) %>%
    left_join(fgeno, by = c("chrom", "pos", "strain2" = "strain")) %>%
    mutate(genotype = as.integer(genotype.x / 2 + genotype.y / 2),
           strain_pair = str_c(strain1, "_", strain2)) %>%
    select(-genotype.x, -genotype.y, -strain1, -strain2)

geno_counts <- geno %>%
    count(chrom, pos, genotype) %>%
    complete(nesting(chrom, pos), genotype = 0:2, fill = list(n = 0)) %>%
    mutate(fraction = n / n_distinct(geno$individual))

fpairs %>%
    # filter(chrom <= 2) %>%
    left_join(geno_counts, by = c("chrom", "pos", "genotype")) %>%
    ggplot(aes(x = pos, y = strain_pair, fill = fraction)) +
    facet_wrap(~ chrom, ncol = 1, scales = "free_x", strip.position = "right") +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +
    theme_minimal() +
    # scale_y_discrete(limits = levels(fgeno$strain)) +
    theme(
        # strip.text = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
    ) +
    ylab("Strain pair")

ggsave("haplotypes/all_chr_pair_sim.png", width = 8, height = 12)

# Same but per individual for one chromosome
tmp <- fpairs %>%
    filter(chrom == 12) %>%
    # filter(as.integer(as.character(pos)) < 1e7) %>%
    rename(founder_geno = genotype) %>%
    left_join(filter(geno, chrom == 12), by = c("chrom", "pos")) %>%
    # filter(individual %in% unique(geno$individual)[1:20]) %>%
    filter(founder_geno == genotype)
ggplot(tmp, aes(x = pos, y = strain_pair)) +
    facet_wrap(~ individual, ncol = 1) +
    geom_tile() +
    # scale_fill_gradient(low = "white", high = "black") +
    theme_minimal() +
    # scale_y_discrete(limits = levels(fgeno$strain)) +
    theme(
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
    ) +
    ylab("Strain pair")

ggsave("haplotypes/pair_match_chr12.png", width = 8, height = 24)

# Shared haplotype regions between founder strains
founder_sim <- strain_pairs %>%
    # filter(chrom <= 2) %>%
    left_join(fgeno, by = c("strain1" = "strain")) %>%
    left_join(fgeno, by = c("chrom", "pos", "strain2" = "strain")) %>%
    mutate(strain_pair = str_c(strain1, "_", strain2),
           geno_match = genotype.x == genotype.y) %>%
    group_by(chrom, strain_pair) %>%
    mutate(roll11_match = slider::slide_dbl(
        geno_match, ~mean(.x, na.rm = TRUE), .before = 5, .after = 5)
        ) %>%
    ungroup() %>%
    select(-genotype.x, -genotype.y)

founder_sim %>%
    # filter(chrom <= 2) %>%
    filter(strain1 != strain2) %>%
    mutate(roll11_match = roll11_match == 1) %>%
    ggplot(aes(x = pos, y = strain_pair, fill = roll11_match)) +
    facet_wrap(~ chrom, ncol = 1, scales = "free_x", strip.position = "right") +
    geom_tile() +
    scale_fill_manual(values = c("white", "black"), na.value = "red") +
    # scale_fill_gradient(low = "white", high = "black") +
    theme_minimal() +
    # scale_y_discrete(limits = levels(fgeno$strain)) +
    theme(
        # strip.text = element_blank(),
        axis.text.y = element_blank(),
        # panel.grid = element_blank(),
        axis.text.x = element_blank(),
    ) +
    ylab("Strain pair") +
    ggtitle("Regions that match between founder pairs")

ggsave("haplotypes/all_chr_founder_sim.png", width = 8, height = 12)
