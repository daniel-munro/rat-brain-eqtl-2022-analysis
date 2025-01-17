library(GenomicRanges)
library(tidyverse)
library(patchwork)

# Include all significant, or all tied top SNPs? The latter would exclude
# independent eQTLs, but maybe that's fine for proportions.
all_top <- tibble(tissue = c("Acbc", "IL", "LHB", "PL", "VoLo")) %>%
    group_by(tissue) %>%
    summarise(
        read_tsv(str_glue("data/tensorqtl/{str_sub(tissue, 1, 1)}QCT.cis_qtl_signif.txt.gz"),
                 col_types = "cc----d---"),
        .groups = "drop"
    ) %>%
    group_by(tissue, phenotype_id) %>%
    filter(pval_nominal == min(pval_nominal)) %>%
    ungroup() %>%
    distinct(tissue, variant_id)

genes <- read_tsv("data/genes.txt", col_types = "c-c---illlll") %>%
    filter(in_expr_Acbc | in_expr_IL | in_expr_LHB | in_expr_PL | in_expr_VoLo)
cis_rng <- with(genes, GRanges(chrom, IRanges(tss - 1e6, tss + 1e6)))

vep_all <- read_tsv("data/vep.txt.gz", comment = "##", na = "-",
                    col_types = "cc-c-cc------c") %>%
    rename(variant_id = `#Uploaded_variation`)
vep_all_rng <- vep_all %>%
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) %>%
    mutate(chrom = str_replace(chrom, "chr", "")) %>%
    with(GRanges(chrom, IRanges(pos, pos)))

vep <- vep_all %>%
    filter(countOverlaps(vep_all_rng, cis_rng) > 0) %>%
    select(variant_id, Consequence) %>%
    separate_rows(Consequence, sep = ",") %>%
    distinct(variant_id, Consequence) %>%
    mutate(
        esnp_Acbc = variant_id %in% all_top$variant_id[all_top$tissue == "Acbc"],
        esnp_IL = variant_id %in% all_top$variant_id[all_top$tissue == "IL"],
        esnp_LHB = variant_id %in% all_top$variant_id[all_top$tissue == "LHB"],
        esnp_PL = variant_id %in% all_top$variant_id[all_top$tissue == "PL"],
        esnp_VoLo = variant_id %in% all_top$variant_id[all_top$tissue == "VoLo"]
    )

# vep %>%
#     group_by(variant_id) %>%
#     summarise(Consequence = str_c(Consequence, collapse = " ")) %>%
#     count(Consequence, sort = TRUE) %>%
#     View()
# # Rank consequences to choose one in the case of multiple annotations:
# cons <- c("intergenic_variant", "intron_variant", )
# # Actually, for proportions I think it makes sense to keep multiple annotations
# # For SNPs.

snp_counts <- vep %>%
    distinct(variant_id, .keep_all = TRUE) %>%
    with(c(Acbc = sum(esnp_Acbc),
           IL = sum(esnp_IL),
           LHB = sum(esnp_LHB),
           PL = sum(esnp_PL),
           VoLo = sum(esnp_VoLo)))

labels <- c(
    "5' UTR" = "5_prime_UTR_variant",
    "3' UTR" = "3_prime_UTR_variant",
    "Missense" = "missense_variant",
    "Synonymous" = "synonymous_variant",
    "Splice region"= "splice_region_variant",
    "Intron" = "intron_variant",
    "Noncoding transcript" = "non_coding_transcript_variant",
    # "Noncoding transcript" = "non_coding_transcript_exon_variant",
    "Intergenic - upstream" = "upstream_gene_variant",
    "Intergenic - downstream" = "downstream_gene_variant",
    "Intergenic - other" = "intergenic_variant"
)

enrich <- vep %>%
    mutate(Consequence = if_else(Consequence == "non_coding_transcript_exon_variant",
                                 "non_coding_transcript_variant",
                                 Consequence)) %>%
    distinct(variant_id, Consequence, .keep_all = TRUE) %>%
    group_by(Consequence) %>%
    summarise(Acbc = sum(esnp_Acbc),
              IL = sum(esnp_IL),
              LHB = sum(esnp_LHB),
              PL = sum(esnp_PL),
              VoLo = sum(esnp_VoLo),
              n_total = n(),
              .groups = "drop") %>%
    pivot_longer(Acbc:VoLo, names_to = "tissue", values_to = "n_esnp") %>%
    mutate(
        #frac_esnp = n_esnp / sum(n_esnp),
        #frac_total = n_total / sum(n_total),
        frac_esnp = n_esnp / snp_counts[tissue],
        frac_total = n_total / n_distinct(vep$variant_id),
        log2_enrich = log2(frac_esnp / frac_total)
    ) %>%
    filter(
        n_total > 100, # Only drops 0.01% of annotations.
        Consequence != "NMD_transcript_variant" # Doesn't describe the variant itself
    ) %>%
    mutate(Consequence = Consequence %>%
               fct_recode(!!!labels) %>%
               fct_relevel(names(labels))) %>%
    arrange(Consequence)

p1 <- enrich %>%
    group_by(Consequence) %>%
    summarise(mean_enr = mean(log2_enrich),
              sd_enr = sd(log2_enrich)) %>%
    mutate(Consequence = fct_rev(Consequence)) %>%
    ggplot(aes(x = Consequence, y = mean_enr, ymin = mean_enr - sd_enr, ymax = mean_enr + sd_enr)) +
    geom_pointrange(fatten = 1) +
    geom_hline(yintercept = 0, lty = 2) +
    coord_flip() +
    theme_minimal() +
    xlab(NULL) +
    ylab(expression(log[2]*" fold enrichment"))

# ggsave("vep/vep_figure_enrich.png", width = 3.5, height = 2)

p2 <- enrich %>%
    group_by(Consequence) %>%
    summarise(mean_frac = mean(frac_esnp),
              sd_frac = sd(frac_esnp)) %>%
    mutate(Consequence = fct_rev(Consequence)) %>%
    ggplot(aes(x = Consequence, y = mean_frac, ymin = mean_frac - sd_frac, ymax = mean_frac + sd_frac)) +
    geom_col(fill = "#aaaaaa") +
    geom_linerange() +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    xlab(NULL) +
    ylab("Proportion\nof variants")

# ggsave("vep/vep_figure_prop.png", width = 1.5, height = 2)

p1 + p2 + plot_layout(widths = c(3, 2))
ggsave("vep/vep_figure.png", width = 5, height = 2.2)

##################################
## Effect size of each category ##
##################################

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciciiidddddidd") %>%
    left_join(select(vep, variant_id, Consequence), by = "variant_id") %>%
    mutate(abs_log_afc = abs(log2_aFC_eQTL))

eqtls %>%
    mutate(Consequence = if_else(Consequence == "non_coding_transcript_exon_variant",
                                 "non_coding_transcript_variant",
                                 Consequence)) %>%
    distinct(variant_id, Consequence, .keep_all = TRUE) %>%
    filter(Consequence %in% labels) %>%
    mutate(Consequence = Consequence %>%
               fct_recode(!!!labels) %>%
               fct_relevel(names(labels)) %>%
               fct_rev()) %>%
    ggplot(aes(x = Consequence, y = abs_log_afc)) +
    geom_boxplot(outlier.size = 0.25) +
    coord_flip(ylim = c(0, 2)) +
    # coord_cartesian(xlim = c(NA, 2)) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank())
