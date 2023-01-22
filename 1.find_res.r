suppressMessages(library(data.table))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

cat("Loading data...")
whole <- readRDS("./data/whole.rds")
cat("Done\n")
cat("Calculating best resolution...")
res_arr <- seq(0.05, 1.0, 0.05)
calc_roc <- function(res) {
    cat(paste0("\tCalculating AUC of res: ", res, "\n"))
    cluster_col <- paste0("RNA_snn_res.", res)
    Idents(whole) <- cluster_col
    markers_roc <- FindAllMarkers(whole, test.use = "roc")
    setDT(markers_roc)
    return(markers_roc[myAUC > 0.7, .N, by = cluster]$N %>% min())
}
res_tbl <- data.table(
    res = res_arr,
    n_min_roc = lapply(res_arr, calc_roc)
)
cat("Done\n")
cat("Saving results...")
fwrite(res_tbl, "results/res_determination.csv")
cat("Done\n")
