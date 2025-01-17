library(tidyverse)
library(patchwork)

expr_log <- read_tsv("../data/expression/ensembl-gene_log2_ComBat_Acbc.bed.gz")
expr_rlog <- read_tsv("../data/expression/ensembl-gene_rlog_ComBat_Acbc.bed.gz")

p1 <- expr_log %>%
    ggplot(aes(x = `00077E67B5`, y = `00077E8336`)) +
    geom_point(size = 0.25) +
    ggtitle("log2")

p2 <- expr_rlog %>%
    ggplot(aes(x = `00077E67B5`, y = `00077E8336`)) +
    geom_point(size = 0.25) +
    ggtitle("rlog")

p1 + p2
ggsave("")
