library(tidyverse)
library(patchwork)

top_lm <- read_tsv("data/gemma/NAcc.lm.assoc.txt.gz",
                   col_types = "c-c--------d") |>
    group_by(gene_id) |>
    sample_n(1) |>
    ungroup()

top_lmm <- read_tsv("data/gemma/NAcc.lmm.assoc.txt.gz",
                    col_types = "c-c----------d--") |>
    group_by(gene_id) |>
    sample_n(1) |>
    ungroup()

pve <- read_tsv("data/gemma/NAcc.grm_pve.pve.txt", col_types = "cdd")

top <- full_join(
    top_lm |> select(gene_id, p_lm = p_wald),
    top_lmm |> select(gene_id, p_lmm = p_wald),
    by = "gene_id"
) |>
    left_join(pve, by = "gene_id")

########################################
## Scatter plot of LMM vs LM p-values ##
########################################

p1_r <- with(top, cor(log10(p_lm), log10(p_lmm))) |> signif(4)
# Colors obtained from viridisLite::viridis(5):
# p1_cols <- c("#440154FF", "#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
# p1_vals <- scales::rescale(c(min(top$pve), c(0, 0.1, 0.4, 0.7, 1) * max(top$pve)))
# p1_cols <- c("#440154FF", "#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF", "#FDE725FF")
# Colors obtained from viridisLite::viridis(6)[1:5] to avoid hard-to-see yellow:
p1_cols <- c("#440154FF", "#440154FF", "#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#7AD151FF")
p1_vals <- scales::rescale(c(min(top$pve), c(0, 0.1, 0.3, 0.5, 0.7, 1) * max(top$pve)))
p1 <- top |>
    mutate(log10p_lm = -log10(p_lm),
           log10p_lmm = -log10(p_lmm)) |>
    slice_sample(prop = 1L, replace = FALSE) |>
    ggplot(aes(x = log10p_lm, y = log10p_lmm, color = pve)) +
    geom_point(size = 0.25) +
    # scale_color_viridis_c(limits = c(0, NA)) +
    # scale_color_viridis_c(limits = c(0, NA), oob = scales::squish) +
    # scale_color_gradient2(low = "red", mid = "gray", high = "blue") +
    scale_color_gradientn(colors = p1_cols, values = p1_vals) +
    annotate("text", x = 2, y = 35, hjust = 0,
             label = str_glue("r = {p1_r}\nN = {nrow(top)}")) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(color = expression(PVE[GRM])) +
    # ggtitle("p-values using LM vs. LMM (NAcc)") +
    xlab(expression(-log[10]*"P, fixed effects model")) +
    ylab(expression(-log[10]*"P, linear mixed model"))
p1

# p1_rinset <- with(top, cor(p_lm, p_lmm)) |> signif(4)
# p1_inset <- top |>
#     slice_sample(prop = 1L, replace = FALSE) |>
#     ggplot(aes(x = p_lm, y = p_lmm, color = pve)) +
#     geom_point(size = 0.02, show.legend = FALSE) +
#     scale_color_gradientn(colors = p1_cols, values = p1_vals) +
#     annotate("text", x = 0.05, y = 1, hjust = 0,
#              label = str_glue("R = {p1_rinset}")) +
#     theme_minimal() +
#     theme(plot.background = element_rect(fill = "white", color = NA),
#           panel.grid = element_blank()) +
#     xlab(NULL) +
#     ylab(NULL)
# p1_inset
# ggsave("GEMMA/gemma_fig_inset.png", width = 2, height = 2)

###########################################
## eGene overlap at different thresholds ##
###########################################

overlap <- tibble(threshold = 10 ^ seq(from = -10, to = -2, length.out = 500)) |>
    group_by(threshold) |>
    summarise(`LM only` = sum(top$p_lm < threshold & top$p_lmm >= threshold),
              `LMM only` = sum(top$p_lm >= threshold & top$p_lmm < threshold),
              both = sum(top$p_lm < threshold & top$p_lmm < threshold),
              .groups = "drop") |>
    pivot_longer(-threshold, names_to = "group", values_to = "eGenes") |>
    mutate(group = fct_relevel(group, "LM only", "both", "LMM only"))

p2 <- overlap |>
    ggplot(aes(x = threshold, y = eGenes, fill = group)) +
    geom_col(position = "fill", width = 0.02) +
    scale_x_log10(expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_fill_manual(values = c("#ff5555", "#aa55aa", "#5555ff")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank()) +
    xlab("P-value threshold") +
    ylab("Proportion of eGenes") +
    labs(fill = "Below\nthreshold in")

p1 / p2 + plot_layout(heights = c(2, 1))
# (p1 + inset_element(p1_inset, 0.52, 0, 1, 0.48)) / p2 + plot_layout(heights = c(2, 1))

ggsave("GEMMA/gemma_fig.png", width = 5, height = 5.5)

#############################
## Stats related to figure ##
#############################

top |>
    summarise(Pearson_p = cor(p_lm, p_lmm),
              Pearson_log_p = cor(log10(p_lm), log10(p_lmm)),
              Spearman_p = cor(p_lm, p_lmm, method = "spearman"))

tibble(threshold = c(1e-9, 1e-6, 1e-3)) |>
    group_by(threshold) |>
    summarise(both = sum(top$p_lm < threshold & top$p_lmm < threshold),
              total = sum(top$p_lm < threshold | top$p_lmm < threshold),
              .groups = "drop") |>
    mutate(fraction = both / total)
