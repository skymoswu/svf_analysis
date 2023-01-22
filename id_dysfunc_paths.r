library(CellChat)
library(tidyverse)
library(circlize)

mod_lst <- c("Mod1", "Mod2", "Mod3")

for (mod_name in mod_lst) {
  cat(sprintf("Processing: %s\n", mod_name))
  # Read in the data
  cc_lst <- readRDS(sprintf("./data/comObjs/%sLst.rds", mod_name))
  cc <- mergeCellChat(cc_lst, add.names = c("LC", "OC", "OD"))
  mod_size <- length(levels(cc@meta$desc))
  oc_obj <- mergeCellChat(cc_lst[c(1, 2)], add.names = c("LC", "OC"))
  od_obj <- mergeCellChat(cc_lst[c(1, 3)], add.names = c("LC", "OD"))
  # Identify DEGs under OC and OD conditions
  oc_net <- identifyOverExpressedGenes(
    oc_obj, group.dataset = "datasets",
    pos.dataset = "OC", features.name = "OC",
    only.pos = FALSE, thresh.pc = 0.1, thresh.p = 0.05
  ) %>% netMappingDEG(., features.name = "OC")
  od_net <- identifyOverExpressedGenes(
    od_obj, group.dataset = "datasets",
    pos.dataset = "OD", features.name = "OD",
    only.pos = FALSE, thresh.pc = 0.1, thresh.p = 0.05
  ) %>% netMappingDEG(., features.name = "OD")
  write.csv(oc_net, sprintf("./results/com/%s_oc.csv", mod_name))
  write.csv(od_net, sprintf("./results/com/%s_od.csv", mod_name))
  # Subset upregulated and downregulated compaths
  oc_net_up <- subsetCommunication(
    oc_obj,
    net = oc_net, datasets = "OC",
    ligand.logFC = 0.2, receptor.logFC = 0.2
  )
  oc_net_down <- subsetCommunication(
    oc_obj,
    net = oc_net, datasets = "LC",
    ligand.logFC = -0.2, receptor.logFC = -0.2
  )
  od_net_up <- subsetCommunication(
    od_obj,
    net = od_net, datasets = "OD",
    ligand.logFC = 0.2, receptor.logFC = 0.2
  )
  od_net_down <- subsetCommunication(
    od_obj,
    net = od_net, datasets = "LC",
    ligand.logFC = -0.2, receptor.logFC = -0.2
  )
  # Plotting for each condition
  # Create path for supplementary figures
  dir.create(sprintf("./plots/com/%s", mod_name))
  png(
    sprintf("./plots/com/%s/oc_up.png", mod_name),
    width = 12, height = 12, res = 300, units = "in"
  )
  netVisual_chord_gene(
    cc_lst[[2]], sources.use = c(1:mod_size),
    targets.use = c(1:mod_size), slot = "net",
    net = oc_net_up,
    title = sprintf("%s: obese control up-regulated", mod_name)
  )
  dev.off()
  png(
    sprintf("./plots/com/%s/oc_down.png", mod_name),
    width = 12, height = 12, res = 300, units = "in"
  )
  netVisual_chord_gene(
    cc_lst[[1]], sources.use = c(1:mod_size),
    targets.use = c(1:mod_size), slot = "net",
    net = oc_net_down,
    title = sprintf("%s: obese control down-regulated", mod_name)
  )
  dev.off()
  png(
    sprintf("./plots/com/%s/od_up.png", mod_name),
    width = 12, height = 12, res = 300, units = "in"
  )
  netVisual_chord_gene(
    cc_lst[[3]], sources.use = c(1:mod_size),
    targets.use = c(1:mod_size), slot = "net",
    net = od_net_up,
    title = sprintf("%s: obese diabetic up-regulated", mod_name)
  )
  dev.off()
  png(
    sprintf("./plots/com/%s/od_down.png", mod_name),
    width = 12, height = 12, res = 300, units = "in"
  )
  netVisual_chord_gene(
    cc_lst[[1]], sources.use = c(1:mod_size),
    targets.use = c(1:mod_size), slot = "net",
    net = od_net_down,
    title = sprintf("%s: obese diabetic down-regulated", mod_name)
  )
  dev.off()
  # Get pathways up/down regulated under both conditions
  up_inter <- intersect(
    oc_net_up$interaction_name,
    od_net_up$interaction_name
  )
  down_inter <- intersect(
    oc_net_down$interaction_name,
    od_net_down$interaction_name
  )
  if (length(up_inter) > 0) {
    png(sprintf("./plots/com/%s_up.png", mod_name),
      width = 12, height = 12, units = "in", res = 300)
    netVisual_chord_gene(
      cc_lst[[3]], sources.use = c(1:mod_size),
      targets.use = c(1:mod_size), slot = "net",
      net = od_net_up[od_net_up$interaction_name %in% up_inter, ],
      title = sprintf("%s: both up-regulated", mod_name)
    )
    dev.off()
  }
  if (length(down_inter) > 0) {
    png(sprintf("./plots/com/%s_down.png", mod_name),
      width = 12, height = 12, units = "in", res = 300)
    netVisual_chord_gene(
      cc_lst[[1]], sources.use = c(1:mod_size),
      targets.use = c(1:mod_size), slot = "net",
      net = od_net_down[od_net_down$interaction_name %in% down_inter, ],
      title = sprintf("%s: both down-regulated", mod_name)
    )
    dev.off()
  }
}
