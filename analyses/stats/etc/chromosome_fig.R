suppressPackageStartupMessages(library(tidyverse))
library(patchwork)

# eqtls <- read_tsv("data/eqtls_indep.txt", col_types = "ccciidddd")
eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciciiidddddidd")

genes <- tibble(tissue = c("Acbc", "IL", "LHB", "PL", "VoLo")) %>%
    group_by(tissue) %>%
    summarise(
        read_tsv(str_glue("data/expression/ensembl-gene_inv-quant_{tissue}.bed.gz"),
                 col_types = cols(`#chr` = "i", start = "i",
                                  gene_id = "c", .default = "-")) %>%
            rename(chrom = `#chr`, pos = start),
        .groups = "drop"
    ) %>%
    distinct(gene_id, chrom, pos)

snps <- read_tsv("data/genotype/P50.rnaseq.88.unpruned.vcf.gz", comment = "##",
                col_types = cols(`#CHROM` = "i", POS = "i", .default = "-")) %>%
    rename(chrom = `#CHROM`, pos = POS)

eqtls %>%
    mutate(chrom = factor(chrom)) %>%
    ggplot(aes(x = chrom, y = pos)) +
    # geom_point()
    geom_violin(scale = "count")

genes %>%
    mutate(chrom = factor(chrom)) %>%
    ggplot(aes(x = chrom, y = pos)) +
    geom_violin(scale = "count")

snps %>%
    mutate(chrom = factor(chrom)) %>%
    ggplot(aes(x = chrom, y = pos)) +
    geom_violin(scale = "count")

genes %>%
    arrange(pos) %>%
    mutate(is_eGene = gene_id %in% eqtls$gene_id) %>%
    ggplot(aes(x = chrom, y = pos, color = is_eGene)) +
    geom_jitter(size = 0.1, alpha = 0.75) +
    theme_minimal()

ggsave("stats/chrom_jitter.png")

genes %>%
    arrange(pos) %>%
    mutate(is_eGene = gene_id %in% eqtls$gene_id,
           seg_end = if_else(is_eGene, 0.4, -0.4)) %>%
    ggplot(aes(x = chrom, xend = chrom + seg_end, y = pos, yend = pos, color = is_eGene)) +
    geom_segment(size = 0.1, alpha = 0.75) +
    theme_minimal()

ggsave("stats/chrom_lines.png")

#######################################
## Count eqtls/genes/snps per Mb bin ##
#######################################

eqtls_per_bin <- eqtls %>%
    mutate(mb = floor(pos / 1e7) * 10) %>%
    count(chrom, mb, name = "n_eqtls")

genes_per_bin <- genes %>%
    mutate(mb = floor(pos / 1e7) * 10) %>%
    count(chrom, mb, name = "n_genes")

snps_per_bin <- snps %>%
    mutate(mb = floor(pos / 1e7) * 10) %>%
    count(chrom, mb, name = "n_snps")

per_mb <- eqtls_per_bin %>%
    full_join(genes_per_bin, by = c("chrom", "mb")) %>%
    full_join(snps_per_bin, by = c("chrom", "mb")) %>%
    replace_na(list(n_eqtls = 0, n_genes = 0, n_snps = 0))

# per_mb %>%
#     pivot_longer(n_eqtls:n_snps, names_to = "feature", values_to = "n") %>%
#     ggplot(aes(x = mb, y = n, color = feature)) +
#     facet_grid(rows = "chrom") +
#     geom_line()

p1 <- per_mb %>%
    ggplot(aes(x = mb, y = n_eqtls)) +
    facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x") +
    geom_line()

p2 <- per_mb %>%
    ggplot(aes(x = mb, y = n_genes)) +
    facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x") +
    geom_line()

p3 <- per_mb %>%
    ggplot(aes(x = mb, y = n_snps)) +
    facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x") +
    geom_line()

p4 <- per_mb %>%
    mutate(eqtls_per_gene = n_eqtls / n_genes) %>%
    ggplot(aes(x = mb, y = eqtls_per_gene)) +
    facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x") +
    geom_line()

p5 <- per_mb %>%
    mutate(eqtls_per_snp = n_eqtls / n_snps) %>%
    ggplot(aes(x = mb, y = eqtls_per_snp)) +
    facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x") +
    geom_line()

p6 <- per_mb %>%
    mutate(eqtls_per_gene_x_snp = n_eqtls / (n_genes * n_snps)) %>%
    ggplot(aes(x = mb, y = eqtls_per_gene_x_snp)) +
    facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x") +
    geom_line()
    # ylim(c(0, 0.0005))

p1 / p2 / p3 / p4 / p5 / p6

ggsave("stats/chromosome_fig.png", width = 12, height = 10)

## Scatter plots for bin counts

q1 <- per_mb %>%
    ggplot(aes(x = n_snps, y = n_genes)) +
    geom_point(size = 0.5, alpha = 0.5)

q2 <- per_mb %>%
    ggplot(aes(x = n_snps, y = n_eqtls)) +
    geom_point(size = 0.5, alpha = 0.5)

q3 <- per_mb %>%
    ggplot(aes(x = n_genes, y = n_eqtls)) +
    geom_point(size = 0.5, alpha = 0.5)

q4 <- per_mb %>%
    mutate(n_snps_x_n_genes = n_snps * n_genes) %>%
    ggplot(aes(x = n_snps_x_n_genes, y = n_eqtls)) +
    geom_point(size = 0.5, alpha = 0.5)

(q1 + q2) / (q3 + q4)
ggsave("stats/chrom_scatter.png")

with(per_mb, cor(n_genes, n_eqtls))
with(per_mb, cor(n_snps, n_eqtls))
