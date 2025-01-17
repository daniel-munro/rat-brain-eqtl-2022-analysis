# Based on https://github.com/PejLab/anevah/blob/master/vignettes/anevah.Rmd

library(tidyverse)

run_anevah <- function(tissue) {
    ref_counts <- read.delim(
        str_glue("data/anevah/anevah_{tissue}_ref_counts.tsv"),
        header = TRUE,
        row.names = 1,
        as.is = TRUE
    )
    alt_counts <- read.delim(
        str_glue("data/anevah/anevah_{tissue}_alt_counts.tsv"),
        header = TRUE,
        row.names = 1,
        as.is = TRUE
    )
    anevah::anevah(ref_counts, alt_counts)
}

for (tissue in c("IL", "LHb", "NAcc", "OFC", "PL")) {
    out <- run_anevah(tissue)
    write_tsv(out, str_glue("analysis/Vg/anevah/Vg.{tissue}.tsv.gz"))
}
