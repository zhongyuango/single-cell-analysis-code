  library(Seurat)
  library(harmony)
  library(DoubletFinder)
  library(SingleR)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)


set.seed(1234)

project_dir <- "."

data_dir   <- file.path(project_dir, "data", "raw")
ref_dir    <- file.path(project_dir, "data", "reference")
result_dir <- file.path(project_dir, "results")

object_dir <- file.path(result_dir, "objects")
table_dir  <- file.path(result_dir, "tables")
figure_dir <- file.path(result_dir, "figures")
qc_dir     <- file.path(figure_dir, "QC")

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,     recursive = TRUE, showWarnings = FALSE)

read_10x_data <- function(input_path,
                          input_type = c("10x_dir", "h5"),
                          matrix_index = NA_integer_) {
  input_type <- match.arg(input_type)

  if (input_type == "10x_dir") {
    counts <- Read10X(data.dir = input_path)
  } else {
    counts <- Read10X_h5(filename = input_path)
  }


add_qc_metrics <- function(object) {
  object[["percent.mt"]] <-
    PercentageFeatureSet(object, pattern = "^MT-")

  object[["percent_ribo"]] <-
    PercentageFeatureSet(object, pattern = "^RP[SL]")

  object[["percent_hb"]] <-
    PercentageFeatureSet(object, pattern = "^HB[^(P)]")

  return(object)
}


plot_qc <- function(object, file_name) {
  p <- VlnPlot(
    object,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "percent_ribo",
      "percent_hb"
    ),
    ncol = 5,
    group.by = "orig.ident",
    pt.size = 0.1
  )

  ggsave(
    filename = file_name,
    plot = p,
    width = 16,
    height = 6,
    dpi = 300
  )

  return(invisible(p))
}


process_one_sample <- function(sample_row) {
  sample_id <- sample_row$sample_id

  message("Processing sample：", sample_id)

   matrix_index <- suppressWarnings(
    as.integer(sample_row$matrix_index)
  )

  counts <- read_10x_data(
    input_path = sample_row$input_path,
    input_type = sample_row$input_type,
    matrix_index = matrix_index
  )

  object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = sample_row$min_cells,
    min.features = sample_row$min_features_create
  )

  object$sample_id <- sample_id
  object$dataset <- sample_row$dataset
  object$condition <- sample_row$condition

  object <- add_qc_metrics(object)

  plot_qc(
    object,
    file.path(qc_dir, paste0(sample_id, "_before_QC.png"))
  )

  object <- subset(
    object,
    subset =
      nFeature_RNA > sample_row$min_features_qc &
      nFeature_RNA < sample_row$max_features_qc &
      percent.mt < sample_row$max_percent_mt
  )

  plot_qc(
    object,
    file.path(qc_dir, paste0(sample_id, "_after_QC.png"))
  )

  object <- NormalizeData(
    object,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  object <- FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = 4000,
    verbose = FALSE
  )

  saveRDS(
    object,
    file.path(object_dir, paste0(sample_id, "_QC.rds"))
  )

  return(object)
}


run_marker_analysis <- function(object,
                                prefix,
                                group_by = "seurat_clusters",
                                min_pct = 0.25,
                                logfc_threshold = 1,
                                top_n = 5) {
  Idents(object) <- object[[group_by, drop = TRUE]]
  DefaultAssay(object) <- "RNA"

  markers <- FindAllMarkers(
    object,
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold
  )

  significant_markers <- markers |>
    filter(p_val_adj < 0.05)

  top_markers <- significant_markers |>
    group_by(cluster) |>
    slice_max(
      order_by = avg_log2FC,
      n = top_n,
      with_ties = FALSE
    ) |>
    ungroup()

  write.csv(
    markers,
    file.path(table_dir, paste0(prefix, "_all_markers.csv")),
    row.names = FALSE
  )

  write.csv(
    significant_markers,
    file.path(table_dir, paste0(prefix, "_significant_markers.csv")),
    row.names = FALSE
  )

  write.csv(
    top_markers,
    file.path(table_dir, paste0(prefix, "_top", top_n, "_markers.csv")),
    row.names = FALSE
  )

  return(list(
    all = markers,
    significant = significant_markers,
    top = top_markers
  ))
}


assign_cluster_labels <- function(object,
                                  cluster_column,
                                  mapping,
                                  output_column) {
  cluster_id <- as.character(
    object@meta.data[[cluster_column]]
  )

  object@meta.data[[output_column]] <-
    unname(mapping[cluster_id])

  object@meta.data[[output_column]][
    is.na(object@meta.data[[output_column]])
  ] <- "Unassigned"

  return(object)
}


sample_info <- read.csv(
  file.path(project_dir, "config", "sample_info.csv"),
  stringsAsFactors = FALSE
)

sample_info_run <- sample_info |>
  filter(include == TRUE)

sample_objects <- map(
  split(sample_info_run, seq_len(nrow(sample_info_run))),
  function(one_row) {
    process_one_sample(as.list(one_row))
  }
)

sample_objects <- compact(sample_objects)
names(sample_objects) <- sample_info_run$sample_id[
  sample_info_run$sample_id %in%
    map_chr(sample_objects, ~ unique(.x$sample_id)[1])
]

qc_summary <- map_dfr(
  sample_objects,
  function(x) {
    data.frame(
      sample_id = unique(x$sample_id)[1],
      dataset = unique(x$dataset)[1],
      condition = unique(x$condition)[1],
      n_cells_after_QC = ncol(x)
    )
  }
)

write.csv(
  qc_summary,
  file.path(table_dir, "QC_summary.csv"),
  row.names = FALSE
)
run_doublet_finder <- function(object,
                               expected_doublet_rate = 0.075,
                               dims = 1:20) {
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, npcs = max(dims), verbose = FALSE)
  object <- FindNeighbors(object, dims = dims, verbose = FALSE)
  object <- FindClusters(
    object,
    resolution = 0.5,
    verbose = FALSE
  )
  object <- RunUMAP(object, dims = dims, verbose = FALSE)

  sweep_result <- paramSweep(
    object,
    PCs = dims,
    sct = FALSE
  )

  sweep_stats <- summarizeSweep(
    sweep_result,
    GT = FALSE
  )

  bcmvn <- find.pK(sweep_stats)

  best_pK <- as.numeric(
    as.character(
      bcmvn$pK[which.max(bcmvn$BCmetric)]
    )
  )

  homotypic_prop <- modelHomotypic(
    object$seurat_clusters
  )

  n_exp <- round(
    expected_doublet_rate * ncol(object)
  )

  n_exp_adjusted <- round(
    n_exp * (1 - homotypic_prop)
  )

  object <- doubletFinder(
    object,
    PCs = dims,
    pN = 0.25,
    pK = best_pK,
    nExp = n_exp_adjusted,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  classification_column <- grep(
    "^DF.classifications",
    colnames(object@meta.data),
    value = TRUE
  )

  singlet_cells <- rownames(object@meta.data)[
    object@meta.data[[classification_column]] ==
      "Singlet"
  ]

  object <- subset(
    object,
    cells = singlet_cells
  )

  return(object)
}

singlet_objects <- imap(
  sample_objects,
  function(object, sample_name) {
    message("DoubletFinder：", sample_name)

    singlet <- run_doublet_finder(object)

    saveRDS(
      singlet,
      file.path(
        object_dir,
        paste0(sample_name, "_singlet.rds")
      )
    )

    return(singlet)
  }
)
scRNA <- merge(
  x = singlet_objects[[1]],
  y = singlet_objects[-1],
  add.cell.ids = names(singlet_objects),
  project = "pan_cancer_scRNA"
)

saveRDS(
  scRNA,
  file.path(object_dir, "scRNA_merged_raw.rds")
)

scRNA <- NormalizeData(
  scRNA,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

scRNA <- FindVariableFeatures(
  scRNA,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = FALSE
)

scRNA <- CellCycleScoring(
  scRNA,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = FALSE
)

scRNA <- ScaleData(
  scRNA,
  vars.to.regress = c(
    "S.Score",
    "G2M.Score",
    "percent_ribo",
    "percent.mt"
  ),
  verbose = FALSE
)

scRNA <- RunPCA(
  scRNA,
  npcs = 50,
  verbose = FALSE
)

scRNA <- RunHarmony(
  object = scRNA,
  group.by.vars = "orig.ident",
  reduction = "pca"
)

scRNA <- RunUMAP(
  scRNA,
  reduction = "harmony",
  dims = 1:50
)

scRNA <- RunTSNE(
  scRNA,
  reduction = "harmony",
  dims = 1:50,
  check_duplicates = FALSE
)

scRNA <- FindNeighbors(
  scRNA,
  reduction = "harmony",
  dims = 1:50
)

scRNA <- FindClusters(
  scRNA,
  resolution = seq(0.1, 1, by = 0.1)
)

scRNA$seurat_clusters <- scRNA$RNA_snn_res.1
Idents(scRNA) <- "seurat_clusters"

saveRDS(
  scRNA,
  file.path(object_dir, "scRNA_harmony.rds")
)

p_sample <- DimPlot(
  scRNA,
  reduction = "umap",
  group.by = "orig.ident",
  raster = FALSE
)

ggsave(
  file.path(figure_dir, "UMAP_by_sample.pdf"),
  p_sample,
  width = 14,
  height = 10
)

p_condition <- DimPlot(
  scRNA,
  reduction = "umap",
  group.by = "condition",
  raster = FALSE
)

ggsave(
  file.path(figure_dir, "UMAP_by_condition.pdf"),
  p_condition,
  width = 9,
  height = 7
)

p_cluster <- DimPlot(
  scRNA,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
)

ggsave(
  file.path(figure_dir, "UMAP_by_cluster.pdf"),
  p_cluster,
  width = 10,
  height = 8
)

celltype_marker_result <- run_marker_analysis(
  object = scRNA,
  prefix = "celltype",
  group_by = "celltype",
  min_pct = 0.25,
  logfc_threshold = 1,
  top_n = 5
)

celltype_top_genes <- unique(
  celltype_marker_result$top$gene
)

scRNA <- ScaleData(
  scRNA,
  features = celltype_top_genes,
  assay = "RNA",
  verbose = FALSE
)

p_celltype_dot <- DotPlot(
  scRNA,
  features = celltype_top_genes,
  group.by = "celltype"
) +
  coord_flip()

ggsave(
  file.path(figure_dir, "DotPlot_celltype_top5_markers.pdf"),
  p_celltype_dot,
  width = 20,
  height = 16
)

p_celltype_heatmap <- DoHeatmap(
  scRNA,
  features = celltype_top_genes,
  group.by = "celltype",
  assay = "RNA",
  label = FALSE
) +
  scale_fill_gradientn(
    colors = c("#2166AC", "white", "#B2182B")
  )

ggsave(
  file.path(figure_dir, "Heatmap_celltype_top5_markers.pdf"),
  p_celltype_heatmap,
  width = 16,
  height = 14
)
major_celltypes <- c(
  "Macrophage/Monocyte",
  "Endothelial_cells",
  "Epithelial_cells",
  "T/NK_cells",
  "Fibroblasts",
  "DC",
  "B_cell"
)

major_celltype_markers <- c(
  "CD3G",
  "CD3E",
  "CD3D",
  "CD14",
  "ITGAM",
  "C1QA",
  "PDGFRA",
  "COL1A2",
  "COL1A1",
  "KRT18",
  "KRT8",
  "EPCAM",
  "KDR",
  "CDH5",
  "PECAM1",
  "ITGAX",
  "FLT3",
  "CD1C",
  "CD79A",
  "MS4A1",
  "CD19"
)

p_major_celltype_marker <- plot_marker_dotplot(
  object = scRNA,
  marker_genes = major_celltype_markers,
  group_by = "celltype",
  group_order = major_celltypes,
  output_file = file.path(
    figure_dir,
    "DotPlot_major_celltype_markers.pdf"
  ),
  width = 15,
  height = 6
)

tnk_celltypes <- c(
  "CTL_CD8",
  "NK-like KLRC2+ T",
  "GammaDelta_T",
  "Proliferating_T",
  "IFN（MX1+）_T",
  "CXCL13+_CD8 T",
  "NK-like FGFBP2+ T",
  "Naive_CD4",
  "Treg",
  "CTL（HSPA1A+）_CD8",
  "Th2",
  "Tem_CD8",
  "CTL（TRAV5+）_CD8",
  "Tfh",
  "Th17/22",
  "Naive_CD8",
  "Eff（TNF+）_CD8"
)

tnk_markers <- c(
  "CXCL13",
  "FOXP3",
  "GATA3",
  "KLRB1",
  "MKI67",
  "IFNG",
  "MX1",
  "TRDC",
  "CXCR5",
  "TRAV5",
  "KLRC2",
  "IL7R",
  "TNF",
  "GZMK",
  "GZMB",
  "HSPA1A",
  "FGFBP2"
)

p_tnk_marker <- plot_marker_dotplot(
  object = tnk_obj,
  marker_genes = tnk_markers,
  group_by = "celltype",
  group_order = tnk_celltypes,
  output_file = file.path(
    figure_dir,
    "DotPlot_TNK_subtype_markers.pdf"
  ),
  width = 15,
  height = 9
)

bcell_celltypes <- c(
  "B_IGHG3+",
  "B_HLA-DRB5+",
  "B_EGR1+",
  "B_IGHA1+",
  "B_TCL1A+",
  "B_IGKV3+",
  "Plasma_cell_MZB1+",
  "B_IGHV4+",
  "B_RRM2+",
  "B_IGKV1+",
  "B_IFIT3+"
)

bcell_markers <- c(
  "SDC1",
  "XBP1",
  "MZB1",
  "BLK",
  "BANK1",
  "IGKV3-15",
  "MS4A1",
  "CD79B",
  "IGHM",
  "IGHV4-31",
  "CD27",
  "IGHG1",
  "IGHG3",
  "IGHA2",
  "JCHAIN",
  "IGHA1",
  "TOP2A",
  "MKI67",
  "RRM2",
  "MX1",
  "ISG15",
  "IFIT3",
  "HLA-DRA",
  "CD74",
  "HLA-DRB5",
  "JUN",
  "FOS",
  "EGR1",
  "IL4R",
  "CCR7",
  "TCL1A"
)

p_bcell_marker <- plot_marker_dotplot(
  object = bcell_obj,
  marker_genes = bcell_markers,
  group_by = "celltype",
  group_order = bcell_celltypes,
  output_file = file.path(
    figure_dir,
    "DotPlot_B_cell_subtype_markers.pdf"
  ),
  width = 20,
  height = 8
)

myeloid_celltypes <- c(
  "Mono_FCN1+",
  "Macro_CCL20+",
  "Macro_C1QA+",
  "Macro_SPP1+",
  "Macro_CXCL10+",
  "Macro_APOE+",
  "Macro_MT1G+"
)

myeloid_markers <- c(
  "HMOX1",
  "MT2A",
  "MT1G",
  "VEGFA",
  "MMP9",
  "SPP1",
  "NFKBIA",
  "IL1B",
  "CCL20",
  "IFIT3",
  "ISG15",
  "CXCL10",
  "LPL",
  "APOC1",
  "APOE",
  "C1QC",
  "C1QB",
  "C1QA",
  "S100A9",
  "S100A8",
  "FCN1"
)

p_myeloid_marker <- plot_marker_dotplot(
  object = myeloid_obj,
  marker_genes = myeloid_markers,
  group_by = "celltype",
  group_order = myeloid_celltypes,
  output_file = file.path(
    figure_dir,
    "DotPlot_Monocyte_Macrophage_subtype_markers.pdf"
  ),
  width = 16,
  height = 7
)

dc_celltypes <- c(
  "CXCL9+_DC",
  "CXCL8+_DC",
  "UBD+_DC",
  "FGL2+_DC"
)

dc_markers <- c(
  "FSCN1",
  "CCL19",
  "UBD",
  "CPVL",
  "CST3",
  "FGL2",
  "IL12B",
  "CXCL10",
  "CXCL9",
  "S100A8",
  "IL1B",
  "CXCL8"
)

p_dc_marker <- plot_marker_dotplot(
  object = dc_obj,
  marker_genes = dc_markers,
  group_by = "celltype",
  group_order = dc_celltypes,
  output_file = file.path(
    figure_dir,
    "DotPlot_DC_subtype_markers.pdf"
  ),
  width = 11,
  height = 5
)

endothelial_celltypes <- c(
  "Endo_ACKR1+",
  "Endo_CCL21",
  "Endo_CXCL12+",
  "Endo_DLK1+",
  "Endo_RGS5+"
)

endothelial_markers <- c(
  "ACTA2",
  "PDGFRB",
  "RGS5",
  "KDR",
  "NOTCH1",
  "DLK1",
  "ANGPT2",
  "VCAM1",
  "CXCL12",
  "LYVE1",
  "PROX1",
  "CCL21",
  "ICAM1",
  "SELE",
  "ACKR1"
)

p_endothelial_marker <- plot_marker_dotplot(
  object = endothelial_obj,
  marker_genes = endothelial_markers,
  group_by = "celltype",
  group_order = endothelial_celltypes,
  output_file = file.path(
    figure_dir,
    "DotPlot_Endothelial_subtype_markers.pdf"
  ),
  width = 13,
  height = 6
)

fibroblast_celltypes <- c(
  "Fibro_STMN1+",
  "Fibro_CD36+",
  "Fibro_MYH11+",
  "Fibro_DPT+",
  "Fibro_MMP1+",
  "Fibro_CXCL14+",
  "Fibro_PI16+"
)

fibroblast_markers <- c(
  "TOP2A",
  "MKI67",
  "STMN1",
  "TAGLN",
  "ACTA2",
  "MYH11",
  "ICOSLG",
  "LTA",
  "MMP3",
  "MMP1",
  "CCL2",
  "CXCL14",
  "LPL",
  "FABP4",
  "CD36",
  "COL3A1",
  "COL1A1",
  "DPT",
  "LUM",
  "DCN",
  "PI16"
)

p_fibroblast_marker <- plot_marker_dotplot(
  object = fibroblast_obj,
  marker_genes = fibroblast_markers,
  group_by = "celltype",
  group_order = fibroblast_celltypes,
  output_file = file.path(
    figure_dir,
    "DotPlot_Fibroblast_subtype_markers.pdf"
  ),
  width = 16,
  height = 7
)


cell_composition <- scRNA@meta.data |>
  count(
    sample_id,
    condition,
    celltype,
    name = "cell_number"
  ) |>
  group_by(sample_id) |>
  mutate(
    proportion = cell_number / sum(cell_number)
  ) |>
  ungroup()

write.csv(
  cell_composition,
  file.path(table_dir, "celltype_composition.csv"),
  row.names = FALSE
)

p_composition <- ggplot(
  cell_composition,
  aes(
    x = sample_id,
    y = proportion,
    fill = celltype
  )
) +
  geom_col(width = 0.85) +
  facet_grid(
    ~ condition,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = NULL,
    y = "Cell proportion",
    fill = "Cell type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave(
  file.path(figure_dir, "Celltype_composition.pdf"),
  p_composition,
  width = 18,
  height = 8
)

saveRDS(
  scRNA,
  file.path(object_dir, "scRNA_final_annotated.rds")
)
