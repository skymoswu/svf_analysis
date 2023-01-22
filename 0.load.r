suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(harmony))

root <- "~/bioresource/SVF_dataset/data"
samples <- list.dirs(root, full.names = FALSE, recursive = FALSE)
samples <- samples[samples != "H06"]
group_remap <- c(
    "H" = "H",
    "ON" = "OC",
    "OT" = "OD"
)

cat("Loading data...\n")

load_data <- function(s) {
    full_path <- paste0(root, "/", s, "/")
    cat(paste0("\tReading sample: ", s, " at path ", full_path, "\n"))
    group_lab <- str_extract(s, "[A-Z]+(?=[0-9])")
    group_lab <- group_remap[group_lab] %>% unname()
    sample_id <- str_extract(s, "(?<=[A-Z])\\d+")
    valid_cells <- read_csv(
        paste0(full_path, "doublets_barcodes.csv"), show_col_types = FALSE
    ) %>% filter(!doubletdetection_doublets)
    mat <- Read10X(data.dir = paste0(full_path, "mtx/"))
    n_cell <- ncol(mat)
    metadata <- data.table(
        sample = rep(paste0(group_lab, sample_id), n_cell),
        group = rep(paste0(group_lab), n_cell)
    )
    rownames(metadata) <- colnames(mat)
    sobj <- CreateSeuratObject(
        counts = mat,
        project = paste0(group_lab, sample_id),
        assay = "RNA",
        min.cells = 5,
        min.features = 500,
        meta.data = metadata
    ) %>% subset(., cells = valid_cells$barcode) %>%
    RenameCells(add.cell.id = paste0(group_lab, sample_id)) %>%
    subset(., subset = HBB == 0 & HBA2 == 0)
    sobj[["mt_percentage"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
    sobj <- subset(
        sobj,
        nFeature_RNA > 500 & mt_percentage < 25
    )
    return(sobj)
}

whole_raw <- lapply(samples, load_data)
cat("Saving raw data...\n")
saveRDS(whole_raw, "data/whole_raw.rds")

cat("Merging...")
whole <- reduce(whole_raw, merge)
rm(whole_raw) # Clean workspace
gc(verbose = FALSE)

whole <- whole %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(features = rownames(.))

res_arr <- seq(0.05, 1.0, 0.05)

whole <- RunPCA(whole, features = VariableFeatures(whole)) %>%
    RunHarmony(., group.by.vars = "sample") %>%
    FindNeighbors(., dims = 1:20, reduction = "harmony") %>%
    FindClusters(., resolution = res_arr) %>%
    RunUMAP(., dims = 1:20)

cat("Saving clustered data...\n")
saveRDS(whole, "data/whole.rds")