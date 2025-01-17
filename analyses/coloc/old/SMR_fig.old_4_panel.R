library(tidyverse)
library(patchwork)

chr_len <- c(282763074, 266435125, 177699992, 184226339, 173707219,
             147991367, 145729302, 133307652, 122095297, 112626471,
             90463843, 52716770, 114033958, 115493446, 111246239,
             90668790, 90843779, 88201929, 62275575, 56205956)
label_locs <- cumsum(c(0, chr_len[1:19])) + chr_len / 2
grid_locs <- cumsum(c(0, chr_len))

genes <- read_tsv("data/genes.txt", col_types = "cc----------")

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "cc---c-----------")

smr <- read_tsv("coloc/SMR.tsv", col_types = "cccdcddd") |>
    filter(tissue == "PL",
           trait == "retrofat")

##########
## eQTL ##
##########

esnps <- read_tsv("data/tensorqtl/PQCT.cis_qtl_signif.txt.gz", col_types = "cc----ddd-") |>
    rename(gene_id = phenotype_id) |>
    filter(gene_id %in% eqtls$gene_id) |>
    left_join(
        smr |>
            select(gene_id, variant_id) |>
            mutate(tested = TRUE),
        by = c("gene_id", "variant_id")
    ) |>
    replace_na(list(tested = FALSE))

df <- esnps |>
    # slice_sample(n = 100000) |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(pval_nominal)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(gpos = pos + cumsum(c(0, chr_len))[chrom],
           tested = if_else(tested, "Tested for\ncolocalization", "Other") |>
               fct_rev()) |>
    arrange(desc(tested), gpos)

p1 <- ggplot(df, aes(x = gpos, y = logp, color = tested, shape = tested)) +
    geom_point(size = 0.5) +
    expand_limits(x = c(0, sum(chr_len))) +
    expand_limits(y = c(0, max(df$logp) * 1.05)) +
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0),
                       minor_breaks = grid_locs) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_color_manual(values = c("#984ea3", "#aaaaaa")) +
    scale_shape_manual(values = c(19, 4)) +
    theme_bw() +
    theme(
        # legend.position = c(0.98, 0.95),
        # legend.justification = c(1, 1),
        # legend.background = element_rect(color = "black", size = 0.2),
        # legend.spacing.y = unit(0, "pt"),
        # legend.key.size = unit(15, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.spacing.x = unit(0, "pt"),
        # legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -5, unit = "pt"),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10](P[eQTL]))) +
    # labs(color = "Tested for\ncolocalization")
    labs(color = NULL, shape = NULL)

# ggsave("coloc/SMR_fig_eQTL.png", width = 7, height = 2)

##########
## GWAS ##
##########

gwas <- read_tsv("data/coloc/adiposity_GWAS/allChr_physiological_retrofat.assoc.txt",
                 col_types = "-c------------d",
                 col_names = c("variant_id", "p_score")) |>
    mutate(tested = variant_id %in% smr$variant_id)

df2 <- gwas |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(p_score)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(gpos = pos + cumsum(c(0, chr_len))[chrom],
           tested = if_else(tested, "Tested for\ncolocalization", "Other") |>
               fct_rev()) |>
    arrange(desc(tested), gpos)

p2 <- ggplot(df2, aes(x = gpos, y = logp, color = tested)) +
    geom_point(size = 0.5) +
    expand_limits(x = c(0, sum(chr_len))) +
    expand_limits(y = c(0, max(df2$logp) * 1.05)) +
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0),
                       minor_breaks = grid_locs) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_color_manual(values = c("black", "#aaaaaa")) +
    theme_bw() +
    theme(
        # legend.position = c(0.98, 0.95),
        # legend.justification = c(1, 1),
        # legend.background = element_rect(color = "black", size = 0.2),
        # legend.spacing.y = unit(0, "pt"),
        # legend.key.size = unit(15, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.spacing.x = unit(0, "pt"),
        # legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -5, unit = "pt"),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10](P[GWAS]))) +
    # labs(color = "Tested for\ncolocalization")
    labs(color = NULL)

# ggsave("coloc/SMR_fig_GWAS.png", width = 7, height = 2)

#########################
## All colocalizations ##
#########################

smr_all <- read_tsv("coloc/SMR.tsv", col_types = "cccdcddd") |>
    group_by(tissue, trait) |>
    mutate(sig = p_SMR < 0.05 / n()) |>
    ungroup()

thresh <- smr_all |>
    group_by(tissue, trait) |>
    summarise(threshold = -log10(0.05 / n()), .groups = "drop") |>
    group_by(tissue) |> # Within tissue they're almost identical
    summarise(threshold = median(threshold)) |>
    ungroup()

df4 <- smr_all |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(p_SMR)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(gpos = pos + cumsum(c(0, chr_len))[chrom]) |>
    arrange(gpos)

p4 <- df4 |>
    filter(sig) |>
    ggplot(aes(x = gpos, y = logp, color = tissue, shape = tissue)) +
    geom_point(aes(color = NULL), data = filter(df4, !sig), size = 0.5, color = "#aaaaaa") +
    geom_point(size = 0.75) +
    geom_hline(aes(yintercept = threshold, color = tissue), data = thresh, size = 0.2, alpha = 0.5) +
    expand_limits(x = c(0, sum(chr_len))) +
    expand_limits(y = c(0, max(df4$logp) * 1.05)) +
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0),
                       minor_breaks = grid_locs) +
    scale_y_continuous(expand = c(0.01, 0)) +
    # scale_color_manual(values = c("#aaaaaa", "black")) +
    scale_color_manual(values = c("#377eb8", "#4daf4a", "#e41a1c", "#ff7f00", "#984ea3")) +
    scale_shape_manual(values = c(1:5)) +
    theme_bw() +
    theme(
        # legend.position = c(0.98, 0.95),
        # legend.justification = c(1, 1),
        # legend.background = element_rect(color = "black", size = 0.2),
        # legend.spacing.y = unit(0, "pt"),
        # legend.key.size = unit(15, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        # legend.spacing.x = unit(0, "pt"),
        # legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -5, unit = "pt"),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10](P[SMR]))) +
    labs(color = "Tissue", shape = "Tissue")

#########################
## SMR for PL/retrofat ##
#########################

df3 <- smr |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(p_SMR)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(gpos = pos + cumsum(c(0, chr_len))[chrom]) |>
    arrange(gpos)

p3 <- ggplot(df3, aes(x = gpos, y = logp)) +
    geom_point(size = 0.5, color = "#984ea3", shape = 5) +
    geom_hline(yintercept = with(thresh, threshold[tissue == "PL"]), size = 0.2, color = "#984ea3") +
    expand_limits(x = c(0, sum(chr_len))) +
    expand_limits(y = c(0, max(df4$logp) * 1.05)) +  # Same limits as all coloc panel
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0),
                       minor_breaks = grid_locs) +
    scale_y_continuous(expand = c(0.01, 0)) +
    # scale_color_manual(values = c("#aaaaaa", "black")) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10](P[SMR])))

# ggsave("coloc/SMR_fig_SMR.png", width = 7, height = 2)

#############
## Combine ##
#############

p1 / p2 / p3 / p4
ggsave("coloc/SMR_fig_A-D.png", width = 9, height = 8)
# * Make sure eQTL isnt subsampled
