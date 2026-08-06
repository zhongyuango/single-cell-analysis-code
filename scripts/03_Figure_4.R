library(Seurat)
library(CellChat)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(NMF)
library(ggalluvial)
library(future)
set.seed(1234)
project_dir <- "."
input_dir  <- file.path(project_dir, "data", "processed")
result_dir <- file.path(project_dir, "results", "CellChat")
figure_dir <- file.path(result_dir, "figures")
table_dir  <- file.path(result_dir, "tables")
object_dir <- file.path(result_dir, "objects")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  T_NK = file.path(input_dir, "T_NK_processed.rds"),
  B_cell = file.path(input_dir, "B_cell_processed.rds"),
  Endothelial = file.path(input_dir, "Endothelial_processed.rds"),
  Fibroblast = file.path(input_dir, "Fibroblast_processed.rds")
)

celltype_col <- "celltype"
sample_col <- "orig.ident"
group_col <- "Type"
species <- "human"

database_categories <- c(
  "Secreted Signaling",
  "ECM-Receptor",
  "Cell-Cell Contact"
)

min_cells <- 10
raw_use <- TRUE
n_workers <- 4

selected_celltypes <- c(
  "B_HLA-DRB5+", "B_IGHA1+", "B_IGHG3+",
  "CTL_CD8", "Tex(CXCL13+)_CD8", "Treg",
  "Endo_ACKR1+", "Endo_CCL21+", "Endo_CXCL12+",
  "Fibro_CXCL13+", "Fibro_DPT+", "Fibro_STMN1+"
)

pathways_to_plot <- c("CD99", "CCL")
outgoing_k <- 4
incoming_k <- 5

read_seurat_object <- function(file_path) {
  object <- readRDS(file_path)
  if (!inherits(object, "Seurat")) {
    stop("The input file is not a Seurat object.：", file_path)
  }
  object
}

standardize_metadata <- function(object,
                                 celltype_column,
                                 sample_column,
                                 group_column) {
  if (!celltype_column %in% colnames(object@meta.data)) {
    stop("The metadata does not contain a cell type column：", celltype_column)
  }

  object$celltype <- as.character(object@meta.data[[celltype_column]])

  if (!group_column %in% colnames(object@meta.data)) {
       sample_parts <- str_split_fixed(object@meta.data[[sample_column]], "_", 2)
    object@meta.data[[group_column]] <- sample_parts[, 2]
  }

  object
}

save_base_pdf <- function(filename, width, height, plot_code) {
  pdf(filename, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  force(plot_code)
}

get_cell_indices <- function(cellchat_object, cell_names) {
  cell_levels <- levels(cellchat_object@idents)
  missing_cells <- setdiff(cell_names, cell_levels)
  matched_cells <- intersect(cell_names, cell_levels)
  match(matched_cells, cell_levels)
}

plot_source_target_bubble <- function(cellchat_object,
                                      source_cells,
                                      target_cells,
                                      output_file,
                                      signaling = NULL,
                                      pair_lr = NULL,
                                      width = 12,
                                      height = 8) {
  source_index <- get_cell_indices(cellchat_object, source_cells)
  target_index <- get_cell_indices(cellchat_object, target_cells)

  save_base_pdf(
    output_file,
    width,
    height,
    {
      print(netVisual_bubble(
        cellchat_object,
        sources.use = source_index,
        targets.use = target_index,
        signaling = signaling,
        pairLR.use = pair_lr,
        remove.isolate = FALSE,
        angle.x = 45
      ))
    }
  )
}

plot_pathway_outputs <- function(cellchat_object,
                                 pathway,
                                 output_directory,
                                 receiver_cells = NULL,
                                 group_size = NULL) {
  receiver_index <- NULL
  if (!is.null(receiver_cells)) {
    receiver_index <- get_cell_indices(cellchat_object, receiver_cells)
  }

  save_base_pdf(file.path(output_directory, paste0(pathway, "_hierarchy.pdf")), 14, 10, {
    netVisual_aggregate(cellchat_object, signaling = pathway,
                        vertex.receiver = receiver_index, layout = "hierarchy")
  })

  save_base_pdf(file.path(output_directory, paste0(pathway, "_circle.pdf")), 9, 8, {
    netVisual_aggregate(cellchat_object, signaling = pathway, layout = "circle")
  })

  save_base_pdf(file.path(output_directory, paste0(pathway, "_chord.pdf")), 10, 8, {
    netVisual_aggregate(cellchat_object, signaling = pathway,
                        layout = "chord", vertex.size = group_size)
  })

  save_base_pdf(file.path(output_directory, paste0(pathway, "_heatmap.pdf")), 8, 7, {
    print(netVisual_heatmap(cellchat_object, signaling = pathway,
                            color.heatmap = "Reds"))
  })

  save_base_pdf(file.path(output_directory, paste0(pathway, "_LR_contribution.pdf")), 7, 5, {
    print(netAnalysis_contribution(cellchat_object, signaling = pathway))
  })

  save_base_pdf(file.path(output_directory, paste0(pathway, "_signaling_role.pdf")), 12, 6, {
    print(netAnalysis_signalingRole_network(
      cellchat_object,
      signaling = pathway,
      width = 12,
      height = 6,
      font.size = 10
    ))
  })

  enriched_lr <- extractEnrichedLR(
    cellchat_object,
    signaling = pathway,
    geneLR.return = FALSE
  )

  if (nrow(enriched_lr) > 0) {
    top_pair <- enriched_lr[1, , drop = FALSE]

    save_base_pdf(file.path(output_directory, paste0(pathway, "_top_LR_circle.pdf")), 10, 8, {
      netVisual_individual(cellchat_object, signaling = pathway,
                           pairLR.use = top_pair, layout = "circle")
    })

    if (!is.null(receiver_index) && length(receiver_index) > 0) {
      save_base_pdf(file.path(output_directory, paste0(pathway, "_top_LR_hierarchy.pdf")), 12, 10, {
        netVisual_individual(cellchat_object, signaling = pathway,
                             pairLR.use = top_pair,
                             vertex.receiver = receiver_index,
                             layout = "hierarchy")
      })
    }
  }

  save_base_pdf(file.path(output_directory, paste0(pathway, "_gene_expression.pdf")), 10, 8, {
    print(plotGeneExpression(cellchat_object, signaling = pathway))
  })

  invisible(enriched_lr)
}

object_list <- lapply(input_files, read_seurat_object)
object_list <- lapply(
  object_list,
  standardize_metadata,
  celltype_column = celltype_col,
  sample_column = sample_col,
  group_column = group_col
)

cellchat_seurat <- merge(
  x = object_list[[1]],
  y = object_list[-1],
  add.cell.ids = names(object_list),
  project = "CellChat_analysis"
)

celltype_count <- as.data.frame(table(cellchat_seurat$celltype))
colnames(celltype_count) <- c("celltype", "cell_number")
write.csv(
  celltype_count,
  file.path(table_dir, "CellChat_input_celltype_counts_before_filtering.csv"),
  row.names = FALSE
)

existing_celltypes <- unique(as.character(cellchat_seurat$celltype))
valid_selected_celltypes <- intersect(selected_celltypes, existing_celltypes)
missing_selected_celltypes <- setdiff(selected_celltypes, existing_celltypes)

cellchat_seurat <- subset(
  cellchat_seurat,
  subset = celltype %in% valid_selected_celltypes
)

Idents(cellchat_seurat) <- "celltype"

saveRDS(
  cellchat_seurat,
  file.path(object_dir, "CellChat_input_Seurat.rds")
)

DefaultAssay(cellchat_seurat) <- "RNA"

data_input <- GetAssayData(cellchat_seurat, assay = "RNA", slot = "data")
metadata_input <- data.frame(
  celltype = cellchat_seurat$celltype,
  row.names = colnames(cellchat_seurat),
  stringsAsFactors = FALSE
)

cellchat <- createCellChat(
  object = data_input,
  meta = metadata_input,
  group.by = "celltype"
)

if (tolower(species) == "human") {
  cellchat_db <- CellChatDB.human
} else if (tolower(species) == "mouse") {
  cellchat_db <- CellChatDB.mouse
} else {
  stop("Species can only be human or mouse.")
}

save_base_pdf(file.path(figure_dir, "CellChat_database_categories.pdf"), 7, 5, {
  showDatabaseCategory(cellchat_db)
})

cellchat_db_use <- subsetDB(cellchat_db, search = database_categories)
cellchat@DB <- cellchat_db_use

future::plan("multisession", workers = n_workers)

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = raw_use)
cellchat <- filterCommunication(cellchat, min.cells = min_cells)

communication_lr <- subsetCommunication(cellchat)
write.csv(
  communication_lr,
  file.path(table_dir, "CellChat_LR_level_communications.csv"),
  row.names = FALSE
)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

communication_pathway <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(
  communication_pathway,
  file.path(table_dir, "CellChat_pathway_level_communications.csv"),
  row.names = FALSE
)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
future::plan("sequential")

group_size <- as.numeric(table(cellchat@idents))

save_base_pdf(file.path(figure_dir, "CellChat_network_interaction_count.pdf"), 12, 10, {
  netVisual_circle(
    cellchat@net$count,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Number of interactions"
  )
})

save_base_pdf(file.path(figure_dir, "CellChat_network_interaction_strength.pdf"), 12, 10, {
  netVisual_circle(
    cellchat@net$weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Interaction strength"
  )
})

network_weight <- cellchat@net$weight
n_groups <- nrow(network_weight)
n_columns <- ceiling(sqrt(n_groups))
n_rows <- ceiling(n_groups / n_columns)

save_base_pdf(file.path(figure_dir, "CellChat_network_each_source.pdf"),
              5 * n_columns, 5 * n_rows, {
  par(mfrow = c(n_rows, n_columns), xpd = TRUE)
  for (i in seq_len(nrow(network_weight))) {
    source_matrix <- matrix(
      0,
      nrow = nrow(network_weight),
      ncol = ncol(network_weight),
      dimnames = dimnames(network_weight)
    )
    source_matrix[i, ] <- network_weight[i, ]
    netVisual_circle(
      source_matrix,
      vertex.weight = group_size,
      weight.scale = TRUE,
      edge.weight.max = max(network_weight),
      title.name = rownames(network_weight)[i]
    )
  }
})

save_base_pdf(file.path(figure_dir, "CellChat_all_interactions_bubble.pdf"), 20, 18, {
  print(netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 45))
})


receiver_cells <- c("CTL_CD8", "Tex(CXCL13+)_CD8", "Treg")
for (pathway in pathways_to_plot) {
  plot_pathway_outputs(
    cellchat,
    pathway,
    figure_dir,
    receiver_cells,
    group_size
  )
}

epithelial_cells <- c("Epithelial_1", "Epithelial_2", "Epithelial_3", "Epithelial_4")
fibroblast_cells <- c("Fibro_CXCL13+", "Fibro_DPT+", "Fibro_STMN1+")
t_cells <- c("CTL_CD8", "Tex(CXCL13+)_CD8", "Treg")

plot_source_target_bubble(
  cellchat,
  epithelial_cells,
  fibroblast_cells,
  file.path(figure_dir, "Epithelial_to_Fibroblast_bubble.pdf")
)

plot_source_target_bubble(
  cellchat,
  fibroblast_cells,
  epithelial_cells,
  file.path(figure_dir, "Fibroblast_to_Epithelial_bubble.pdf")
)

plot_source_target_bubble(
  cellchat,
  t_cells,
  fibroblast_cells,
  file.path(figure_dir, "Tcell_to_Fibroblast_bubble.pdf")
)

plot_source_target_bubble(
  cellchat,
  fibroblast_cells,
  t_cells,
  file.path(figure_dir, "Fibroblast_to_Tcell_bubble.pdf")
)

plot_source_target_bubble(
  cellchat,
  fibroblast_cells,
  t_cells,
  file.path(figure_dir, "Fibroblast_to_Tcell_CCL_CXCL_bubble.pdf"),
  signaling = c("CCL", "CXCL")
)

p_role_scatter <- netAnalysis_signalingRole_scatter(
  cellchat,
  label.size = 4,
  font.size = 10
)

ggsave(
  file.path(figure_dir, "CellChat_signaling_role_scatter.pdf"),
  plot = p_role_scatter,
  width = 9,
  height = 7
)

save_base_pdf(file.path(figure_dir, "CellChat_outgoing_signaling_role_heatmap.pdf"), 9, 12, {
  print(netAnalysis_signalingRole_heatmap(
    cellchat,
    pattern = "outgoing",
    width = 10,
    height = 18,
    title = "Outgoing signaling",
    font.size = 10,
    font.size.title = 12
  ))
})

save_base_pdf(file.path(figure_dir, "CellChat_incoming_signaling_role_heatmap.pdf"), 9, 12, {
  print(netAnalysis_signalingRole_heatmap(
    cellchat,
    pattern = "incoming",
    width = 10,
    height = 18,
    title = "Incoming signaling",
    font.size = 10,
    font.size.title = 12
  ))
})

save_base_pdf(file.path(figure_dir, "CellChat_outgoing_selectK.pdf"), 8, 6, {
  selectK(cellchat, pattern = "outgoing")
})

cellchat <- identifyCommunicationPatterns(
  cellchat,
  pattern = "outgoing",
  k = outgoing_k,
  width = 8,
  height = 14
)

outgoing_pattern <- cellchat@netP$pattern$outgoing
write.csv(outgoing_pattern$pattern$cell,
          file.path(table_dir, "CellChat_outgoing_pattern_cells.csv"))
write.csv(outgoing_pattern$pattern$signaling,
          file.path(table_dir, "CellChat_outgoing_pattern_signaling.csv"))

save_base_pdf(file.path(figure_dir, "CellChat_outgoing_river.pdf"), 12, 12, {
  print(netAnalysis_river(
    cellchat,
    pattern = "outgoing",
    color.use.signaling = "C7",
    font.size = 3
  ))
})

save_base_pdf(file.path(figure_dir, "CellChat_outgoing_dotplot.pdf"), 14, 12, {
  print(netAnalysis_dot(cellchat, pattern = "outgoing"))
})

save_base_pdf(file.path(figure_dir, "CellChat_incoming_selectK.pdf"), 8, 6, {
  selectK(cellchat, pattern = "incoming")
})

cellchat <- identifyCommunicationPatterns(
  cellchat,
  pattern = "incoming",
  k = incoming_k,
  width = 8,
  height = 14
)

incoming_pattern <- cellchat@netP$pattern$incoming
write.csv(incoming_pattern$pattern$cell,
          file.path(table_dir, "CellChat_incoming_pattern_cells.csv"))
write.csv(incoming_pattern$pattern$signaling,
          file.path(table_dir, "CellChat_incoming_pattern_signaling.csv"))

save_base_pdf(file.path(figure_dir, "CellChat_incoming_river.pdf"), 12, 12, {
  print(netAnalysis_river(
    cellchat,
    pattern = "incoming",
    color.use.signaling = "C7",
    font.size = 3
  ))
})

save_base_pdf(file.path(figure_dir, "CellChat_incoming_dotplot.pdf"), 14, 12, {
  print(netAnalysis_dot(cellchat, pattern = "incoming"))
})

saveRDS(cellchat, file.path(object_dir, "CellChat.rds"))
saveRDS(cellchat_db_use, file.path(object_dir, "CellChat_database.rds"))

writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "CellChat_sessionInfo.txt")
)
