library(VariantAnnotation)
library(tidyverse)

genes <- tibble(tissue = c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("data/expression/ensembl-gene_log2_{tissue}.bed.gz"),
                 col_types = cols(gene_id = "c", .default = "-")),
        .groups = "drop"
    ) |>
    distinct(gene_id) |>
    left_join(
        read_tsv("data/genes.txt", col_types = "c-c---i-----"),
        by = "gene_id"
    )

###########################
## Count cis-window SNPs ##
###########################

# # This doesn't group identical-genotype SNPs, so I'm using the method below it instead.
# count_cis_snps <- function(chromo, tss) {
#     geno |>
#         filter(chrom == chromo,
#                pos > tss - 1e6,
#                pos < tss + 1e6) |>
#         nrow()
# }
# geno <- info(readVcf("data/genotype/P50.rnaseq.88.unpruned.vcf.gz")) |>
#     as_tibble(rownames = "variant_id") |>
#     mutate(AC = unlist(AC)) |>
#     filter(AC > 0, AC < 88 * 2) |>
#     separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
#     mutate(chrom = str_replace(chrom, "chr", ""))

# Based on locuszoom_fig.R
load_geno <- function(filename, chrom, start, end) {
    rng <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))
    gt <- VariantAnnotation::readGT(filename, param = rng)
    # gt <- VariantAnnotation::readGT(filename)
    geno <- apply(gt, 2, function(x) c("0|0" = 0L, "0|1" = 1L, "1|0" = 1L, "1|1" = 2L)[x])
    geno
}

count_unique_cis_snps <- function(chrom, tss) {
    # print(c(chrom, tss))
    filename <- "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
    geno <- load_geno(filename, chrom, tss - 1e6L, tss + 1e6L)
    if (!is.matrix(geno)) {
        if (length(geno) > 0) {
            geno <- matrix(geno, nrow = 1)
        } else {
            return(0)
        }
    }
    geno |>
        unique() |> # Count sets of identical-genotype SNPs
        apply(1, function(x) length(unique(x)) > 1) |> # Ignore monomorphic (or all het)
        sum()
}

cis_snps <- genes |>
    # dplyr::slice(1:1000) |>
    rowwise() |>
    mutate(n_cis_snps = count_unique_cis_snps(chrom, tss)) |>
    ungroup()

cis_snps |>
    summarise(n_0_snps = sum(n_cis_snps == 0),
              n_10_snps = sum(n_cis_snps < 10),
              n_100_snps = sum(n_cis_snps < 100),
              n_genes = n(),
              frac_0_snps = mean(n_cis_snps == 0),
              frac_10_snps = mean(n_cis_snps < 10),
              frac_100_snps = mean(n_cis_snps < 100))

# n_0_snps n_10_snps n_100_snps n_genes frac_0_snps frac_10_snps frac_100_snps
#       48       596      12569   17487     0.00274       0.0341         0.719

###############################
## Other genes in cis-window ##
###############################

# genes_in_window <- genes |>
#     rowwise() |>
#     mutate(n_in_window = sum(
#         genes$chrom == chrom & genes$tss > tss - 1e6 & genes$tss < tss + 1e6
#     ) - 1L) |>
#     ungroup()
# 
# summary(genes_in_window$n_in_window)

# Count only genes with TSS in high/perfect LD with a gene's TSS

load_geno_all <- function(filename) {
    gt <- VariantAnnotation::readGT(filename)
    geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)[x])
    rownames(geno) <- rownames(gt)
    # geno <- geno[apply(geno, 1, function(x) sd(x) > 0), ]
    rsums <- rowSums(geno)
    geno <- geno[rsums > 0 & rsums < ncol(geno) * 2, ] # Remove monomorphic
    geno
}

nearest_snp <- function(tss_chrom, tss_pos, snps) {
    nearest <- snps |>
        # filter(chrom == tss_chrom) |>
        mutate(dist = abs(pos - tss_pos)) |>
        arrange(dist) |>
        slice(1)
    stopifnot(nearest$chrom == tss_chrom)
    nearest$variant_id
}

geno <- load_geno_all("data/genotype/P50.rnaseq.88.unpruned.vcf.gz")
snps <- tibble(variant_id = rownames(geno)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE, remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", ""))
snps <- split(snps, snps$chrom)
genes2 <- genes |>
    # slice_sample(n = 1000) |>
    rowwise() |>
    mutate(nearest = nearest_snp(chrom, tss, snps[[chrom]])) |>
    ungroup()
# geno2 <- genes2 |>
#     split(genes2$chrom) |>
#     map(function(df) geno[df$nearest, ])
genes_in_high_ld <- genes2 |>
    group_by(chrom) |>
    summarise({
        geno2 <- geno[nearest, ]
        rownames(geno2) <- gene_id # nearest SNPs are not unique
        r2 <- cor(t(geno2)) ^ 2
        r2 |>
            as_tibble(rownames = "gene_id") |>
            pivot_longer(-gene_id, names_to = "gene_id2", values_to = "r2")
    }, .groups = "drop") |>
    filter(gene_id != gene_id2) |>
    group_by(gene_id) |>
    summarise(n_genes_ld_100 = sum(r2 == 1),
              n_genes_ld_99 = sum(r2 > 0.99),
              n_genes_ld_95 = sum(r2 > 0.95),
              n_genes_ld_90 = sum(r2 > 0.9),
              .groups = "drop")
genes_in_high_ld |>
    summarise(across(-gene_id, ~ sum(. > 0)),
              total = n())

# n_genes_ld_100 n_genes_ld_99 n_genes_ld_95 n_genes_ld_90 total
#          10187         12752         14191         15085 17487
