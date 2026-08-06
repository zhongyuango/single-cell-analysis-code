set.seed(1234)

project_dir <- "."

data_dir <- file.path(project_dir, "data")
reference_dir <- file.path(data_dir, "reference")
processed_dir <- file.path(data_dir, "processed")

result_dir <- file.path(project_dir, "results")
copykat_dir <- file.path(result_dir, "copyKAT")
gsva_dir <- file.path(result_dir, "GSVA")
hdwgcna_dir <- file.path(result_dir, "hdWGCNA")
survival_dir <- file.path(result_dir, "survival")

dir.create(copykat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsva_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hdwgcna_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(survival_dir, recursive = TRUE, showWarnings = FALSE)

library(copykat)
library(Seurat)
library(ggplot2)
library(dplyr)

copykat_input_file <- file.path(
  processed_dir,
  "Epithelial_processed.rds"
)

copykat_sample_name <- "Epithelial_cells"

copykat_celltype_col <- "celltype"

normal_cell_labels <- c(
  "AT1",
  "AT2",
  "Ciliated"
)

copykat_ks_cut <- 0.05
copykat_n_gene_chr <- 5
copykat_window_size <- 25
copykat_distance <- "euclidean"
copykat_genome <- "hg20"
copykat_n_cores <- 1


seurat_obj <- readRDS(copykat_input_file)


DefaultAssay(seurat_obj) <- "RNA"

expr_matrix <- tryCatch(
  {
    GetAssayData(
      object = seurat_obj,
      assay = "RNA",
      layer = "counts"
    )
  },
  error = function(e) {
    GetAssayData(
      object = seurat_obj,
      assay = "RNA",
      slot = "counts"
    )
  }
)

normal_cell_names <- NULL

if (copykat_celltype_col %in% colnames(seurat_obj@meta.data)) {
  normal_cell_names <- colnames(seurat_obj)[
    seurat_obj@meta.data[[copykat_celltype_col]] %in% normal_cell_labels
  ]

copykat_result <- copykat(
  rawmat = expr_matrix,
  id.type = "S",
  ngene.chr = copykat_n_gene_chr,
  win.size = copykat_window_size,
  KS.cut = copykat_ks_cut,
  sam.name = file.path(
    copykat_dir,
    copykat_sample_name
  ),
  distance = copykat_distance,
  norm.cell.names = normal_cell_names,
  output.seg = TRUE,
  plot.genes = TRUE,
  genome = copykat_genome,
  n.cores = copykat_n_cores
)

prediction_df <- as.data.frame(copykat_result$prediction)

if ("cell.names" %in% colnames(prediction_df)) {
  prediction_cells <- prediction_df$cell.names
} else {
  prediction_cells <- rownames(prediction_df)
}

copykat_status <- prediction_df$copykat.pred
names(copykat_status) <- prediction_cells

seurat_obj$copykat_status <- copykat_status[colnames(seurat_obj)]

tumor_cells <- names(copykat_status)[
  copykat_status == "aneuploid"
]

normal_cells <- names(copykat_status)[
  copykat_status == "diploid"
]

not_defined_cells <- names(copykat_status)[
  is.na(copykat_status) |
    copykat_status == "not.defined"
]

copykat_summary <- data.frame(
  category = c(
    "Total predicted cells",
    "Aneuploid cells",
    "Diploid cells",
    "Not defined cells"
  ),
  n_cells = c(
    length(copykat_status),
    length(tumor_cells),
    length(normal_cells),
    length(not_defined_cells)
  )
)

print(copykat_summary)
print(table(seurat_obj$copykat_status, useNA = "ifany"))

write.csv(
  prediction_df,
  file = file.path(
    copykat_dir,
    "copyKAT_predictions.csv"
  ),
  row.names = FALSE
)

write.csv(
  copykat_summary,
  file = file.path(
    copykat_dir,
    "copyKAT_summary.csv"
  ),
  row.names = FALSE
)

if ("umap" %in% names(seurat_obj@reductions)) {
  p_copykat <- DimPlot(
    object = seurat_obj,
    reduction = "umap",
    group.by = "copykat_status",
    raster = FALSE
  ) +
    ggtitle("copyKAT classification") +
    theme_classic()

  ggsave(
    filename = file.path(
      copykat_dir,
      "copyKAT_classification_UMAP.pdf"
    ),
    plot = p_copykat,
    width = 8,
    height = 6
  )
}

saveRDS(
  object = copykat_result,
  file = file.path(
    copykat_dir,
    "copyKAT_result.rds"
  )
)

saveRDS(
  object = seurat_obj,
  file = file.path(
    copykat_dir,
    "epithelial_with_copyKAT.rds"
  )
)

tumor_cells <- intersect(
  tumor_cells,
  colnames(seurat_obj)
)

if (length(tumor_cells) > 0) {
  tumor_object <- subset(
    x = seurat_obj,
    cells = tumor_cells
  )

  saveRDS(
    object = tumor_object,
    file = file.path(
      copykat_dir,
      "copyKAT_aneuploid_cells.rds"
    )
  )
}

library(Seurat)
library(GSVA)
library(GSEABase)
library(limma)
library(pheatmap)
library(stringr)
library(dplyr)

gsva_input_file <- file.path(
  copykat_dir,
  "epithelial_with_copyKAT.rds"
)

gmt_file <- file.path(
  reference_dir,
  "c2.cp.reactome.v2024.1.Hs.symbols.gmt"
)

selected_pathway_file <- file.path(
  reference_dir,
  "all_pathway.csv"
)

gsva_group_col <- "Type"

group_1 <- "Immunotherapy"
group_2 <- "Naive"

max_cells_for_gsva <- 5000

gsva_min_size <- 20
gsva_max_size <- 1000

gsva_fdr_cutoff <- 0.05
gsva_logfc_cutoff <- 0

gsva_object <- readRDS(gsva_input_file)

gene_sets <- getGmt(gmt_file)

if (!is.null(max_cells_for_gsva) &&
    ncol(gsva_object) > max_cells_for_gsva) {
  sampled_cells <- sample(
    x = colnames(gsva_object),
    size = max_cells_for_gsva
  )

  gsva_object <- subset(
    x = gsva_object,
    cells = sampled_cells
  )
}

expression_matrix <- tryCatch(
  {
    GetAssayData(
      object = gsva_object,
      assay = "RNA",
      layer = "data"
    )
  },
  error = function(e) {
    GetAssayData(
      object = gsva_object,
      assay = "RNA",
      slot = "data"
    )
  }
)

gsva_scores <- tryCatch(
  {
    gsva_parameter <- gsvaParam(
      exprData = as.matrix(expression_matrix),
      geneSets = gene_sets,
      minSize = gsva_min_size,
      maxSize = gsva_max_size,
      kcdf = "Gaussian"
    )

    gsva(
      param = gsva_parameter,
      verbose = TRUE
    )
  },
  error = function(e) {
    gsva(
      expr = as.matrix(expression_matrix),
      gset.idx.list = gene_sets,
      min.sz = gsva_min_size,
      max.sz = gsva_max_size,
      kcdf = "Gaussian",
      verbose = TRUE
    )
  }
)

write.csv(
  x = gsva_scores,
  file = file.path(
    gsva_dir,
    "GSVA_cell_level_scores.csv"
  )
)

common_cells <- intersect(
  colnames(gsva_scores),
  rownames(gsva_object@meta.data)
)

gsva_scores <- gsva_scores[
  ,
  common_cells,
  drop = FALSE
]

groups <- factor(
  gsva_object@meta.data[
    common_cells,
    gsva_group_col
  ]
)

design <- model.matrix(~ 0 + groups)

colnames(design) <- levels(groups)

contrast_string <- paste0(
  "`",
  group_1,
  "`-`",
  group_2,
  "`"
)

contrast_matrix <- makeContrasts(
  contrasts = contrast_string,
  levels = design
)

fit <- lmFit(
  object = gsva_scores,
  design = design
)

fit <- contrasts.fit(
  fit = fit,
  contrasts = contrast_matrix
)

fit <- eBayes(fit)

gsva_differential <- topTable(
  fit = fit,
  number = Inf,
  adjust.method = "BH"
)

gsva_differential$pathway <- rownames(gsva_differential)

significant_pathways <- gsva_differential |>
  filter(
    adj.P.Val < gsva_fdr_cutoff,
    abs(logFC) > gsva_logfc_cutoff
  ) |>
  arrange(adj.P.Val)

write.csv(
  x = gsva_differential,
  file = file.path(
    gsva_dir,
    paste0(
      "GSVA_",
      group_1,
      "_vs_",
      group_2,
      "_all_pathways.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  x = significant_pathways,
  file = file.path(
    gsva_dir,
    paste0(
      "GSVA_",
      group_1,
      "_vs_",
      group_2,
      "_significant_pathways.csv"
    )
  ),
  row.names = FALSE
)

group_levels <- levels(groups)

average_gsva <- sapply(
  X = group_levels,
  FUN = function(group_name) {
    rowMeans(
      gsva_scores[
        ,
        groups == group_name,
        drop = FALSE
      ]
    )
  }
)

average_gsva <- as.matrix(average_gsva)

write.csv(
  x = average_gsva,
  file = file.path(
    gsva_dir,
    "GSVA_average_scores_by_group.csv"
  )
)

heatmap_matrix <- average_gsva

normalize_pathway_name <- function(x) {
  x |>
    toupper() |>
    gsub(
      pattern = "[- ]",
      replacement = "_"
    )
}

if (file.exists(selected_pathway_file)) {
  selected_pathways <- read.csv(
    file = selected_pathway_file,
    row.names = 1,
    check.names = FALSE
  )

  selected_names <- normalize_pathway_name(
    rownames(selected_pathways)
  )

  significant_names <- normalize_pathway_name(
    significant_pathways$pathway
  )

  matrix_names <- normalize_pathway_name(
    rownames(average_gsva)
  )

  selected_significant <- intersect(
    selected_names,
    significant_names
  )

  keep_index <- matrix_names %in% selected_significant

  heatmap_matrix <- average_gsva[
    keep_index,
    ,
    drop = FALSE
  ]
}

rownames(heatmap_matrix) <- label_mapping$display_label

  write.csv(
    x = label_mapping,
    file = file.path(
      gsva_dir,
      "GSVA_heatmap_pathway_label_mapping.csv"
    ),
    row.names = FALSE
  )

  heatmap_height <- max(
    6,
    0.35 * nrow(heatmap_matrix) + 2
  )

  pdf(
    file = file.path(
      gsva_dir,
      "GSVA_pathway_heatmap.pdf"
    ),
    width = 9,
    height = heatmap_height,
    useDingbats = FALSE
  )

  pheatmap(
    mat = heatmap_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    color = colorRampPalette(
      c(
        "#2166AC",
        "#478ABF",
        "#90C0DC",
        "white",
        "#EF8C65",
        "#CF4F45",
        "#B2182B"
      )
    )(100),
    border_color = NA,
    fontsize_row = 9,
    fontsize_col = 10,
    main = "GSVA analysis"
  )

  dev.off()

  write.csv(
    x = heatmap_matrix,
    file = file.path(
      gsva_dir,
      "GSVA_pathway_heatmap_matrix.csv"
    )
  )
}

library(hdWGCNA)
library(WGCNA)
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(cowplot)
library(msigdbr)
library(UCell)

hdwgcna_input_file <- file.path(
  copykat_dir,
  "copyKAT_aneuploid_cells.rds"
)

hdwgcna_name <- "Tumor"

hdwgcna_group_col <- "Type"

metacell_group_cols <- c(
  "Type"
)

metacell_k <- 25
metacell_max_shared <- 20
metacell_reduction <- "umap"

soft_power <- 5
network_type <- "signed"
tom_name <- "TOM"

n_threads <- 10
n_hub_genes <- 25

module_colors <- c(
  "#E87B89",
  "#1F6E9C",
  "#2B9B81"
)
clinical_info_file <- file.path(
  reference_dir,
  "clinical_information.csv"
)

hdwgcna_object <- readRDS(hdwgcna_input_file)
 metadata <- hdwgcna_object@meta.data |>
    rownames_to_column("cell_id") |>
    left_join(
      clinical_info,
      by = "orig.ident"
    ) |>
    column_to_rownames("cell_id")

  metadata <- metadata[
    colnames(hdwgcna_object),
    ,
    drop = FALSE
  ]

  hdwgcna_object@meta.data <- metadata
}

theme_set(theme_cowplot())

enableWGCNAThreads(
  nThreads = n_threads
)

hdwgcna_object <- SetupForWGCNA(
  hdwgcna_object,
  wgcna_name = hdwgcna_name,
  gene_select = "custom",
  gene_list = valid_genes
)

hdwgcna_object <- MetacellsByGroups(
  seurat_obj = hdwgcna_object,
  group.by = metacell_group_cols,
  k = metacell_k,
  max_shared = metacell_max_shared,
  reduction = metacell_reduction,
  ident.group = hdwgcna_group_col
)

hdwgcna_object <- NormalizeMetacells(
  hdwgcna_object
)

hdwgcna_object <- SetDatExpr(
  hdwgcna_object,
  assay = "RNA",
  slot = "data"
)

hdwgcna_object <- TestSoftPowers(
  hdwgcna_object,
  networkType = network_type
)

soft_power_plots <- PlotSoftPowers(
  hdwgcna_object
)

p_soft_power <- wrap_plots(
  soft_power_plots,
  ncol = 2
)

ggsave(
  filename = file.path(
    hdwgcna_dir,
    "hdWGCNA_soft_power_selection.pdf"
  ),
  plot = p_soft_power,
  width = 8,
  height = 8
)

hdwgcna_object <- ConstructNetwork(
  hdwgcna_object,
  soft_power = soft_power,
  tom_name = tom_name,
  setDatExpr = FALSE,
  overwrite_tom = TRUE
)

pdf(
  file = file.path(
    hdwgcna_dir,
    "hdWGCNA_dendrogram_initial.pdf"
  ),
  width = 10,
  height = 6,
  useDingbats = FALSE
)

PlotDendrogram(
  hdwgcna_object,
  main = "hdWGCNA dendrogram"
)

dev.off()

hdwgcna_object <- ScaleData(
  object = hdwgcna_object,
  verbose = FALSE
)

hdwgcna_object <- ModuleEigengenes(
  seurat_obj = hdwgcna_object,
  group.by.vars = hdwgcna_group_col
)

harmonized_mes <- GetMEs(
  hdwgcna_object,
  harmonized = TRUE
)

original_mes <- GetMEs(
  hdwgcna_object,
  harmonized = FALSE
)

hdwgcna_object <- ModuleConnectivity(
  hdwgcna_object
)

modules <- GetModules(hdwgcna_object)

write.csv(
  x = modules,
  file = file.path(
    hdwgcna_dir,
    "hdWGCNA_modules.csv"
  ),
  row.names = FALSE
)

hub_genes <- GetHubGenes(
  hdwgcna_object,
  n_hubs = n_hub_genes
)

write.csv(
  x = hub_genes,
  file = file.path(
    hdwgcna_dir,
    "hdWGCNA_hub_genes.csv"
  ),
  row.names = FALSE
)

hdwgcna_object <- ModuleExprScore(
  hdwgcna_object,
  n_genes = n_hub_genes,
  method = "UCell"
)

p_kme <- PlotKMEs(
  hdwgcna_object,
  ncol = 3
)

ggsave(
  filename = file.path(
    hdwgcna_dir,
    "hdWGCNA_module_kME.pdf"
  ),
  plot = p_kme,
  width = 8,
  height = 5
)

module_feature_plots <- ModuleFeaturePlot(
  hdwgcna_object,
  reduction = "umap",
  features = "hMEs",
  order = TRUE,
  raster = TRUE
)

p_module_features <- wrap_plots(
  module_feature_plots,
  ncol = 3
)

ggsave(
  filename = file.path(
    hdwgcna_dir,
    "hdWGCNA_hME_UMAP.pdf"
  ),
  plot = p_module_features,
  width = 10,
  height = 6
)

harmonized_mes <- GetMEs(
  hdwgcna_object,
  harmonized = TRUE
)

common_meta_cells <- intersect(
  rownames(hdwgcna_object@meta.data),
  rownames(harmonized_mes)
)

hdwgcna_object@meta.data[
  common_meta_cells,
  colnames(harmonized_mes)
] <- harmonized_mes[
  common_meta_cells,
  ,
  drop = FALSE
]

module_network_dir <- file.path(
  hdwgcna_dir,
  "module_networks"
)

dir.create(
  module_network_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ModuleNetworkPlot(
  hdwgcna_object,
  outdir = module_network_dir
)

if (requireNamespace("enrichR", quietly = TRUE)) {
  enrichr_databases <- c(
    "GO_Biological_Process_2023"
  )

  hdwgcna_object <- RunEnrichr(
    hdwgcna_object,
    dbs = enrichr_databases,
    max_genes = 100
  )

  enrich_df <- GetEnrichrTable(
    hdwgcna_object
  )

  write.csv(
    x = enrich_df,
    file = file.path(
      hdwgcna_dir,
      "hdWGCNA_GO_enrichment.csv"
    ),
    row.names = FALSE
  )

  pdf(
    file = file.path(
      hdwgcna_dir,
      "hdWGCNA_module_GO_dotplot.pdf"
    ),
    width = 8,
    height = 8,
    useDingbats = FALSE
  )

  EnrichrDotPlot(
    hdwgcna_object,
    mods = "all",
    database = "GO_Biological_Process_2023",
    n_terms = 3
  )

  dev.off()
} 

saveRDS(
  object = hdwgcna_object,
  file = file.path(
    hdwgcna_dir,
    "hdWGCNA_object.rds"
  )
)

library(survival)
library(survminer)
library(ggplot2)

survival_input_file <- file.path(
  processed_dir,
  "survival_data.csv"
)

survival_gene <- "ZNF385A"

survival_time_col <- "OS_time"
survival_status_col <- "OS_status"

survival_time_unit <- "days"
survival_output_unit <- "years"

minimum_group_proportion <- 0.10

maximum_followup_years <- 2.5

survival_data <- read.csv(
  file = survival_input_file,
  check.names = FALSE
)

required_columns <- c(
  survival_time_col,
  survival_status_col,
  survival_gene
)

missing_columns <- setdiff(
  required_columns,
  colnames(survival_data)
)

survival_data <- survival_data[
  complete.cases(
    survival_data[
      ,
      required_columns
    ]
  ),
  ,
  drop = FALSE
]

cut_result <- surv_cutpoint(
  data = survival_data,
  time = survival_time_col,
  event = survival_status_col,
  variables = survival_gene,
  minprop = minimum_group_proportion
)

categorized_data <- surv_categorize(
  cut_result
)

group_counts <- as.data.frame(
  table(
    categorized_data[[survival_gene]]
  )
)

colnames(group_counts) <- c(
  "group",
  "n_samples"
)

write.csv(
  x = group_counts,
  file = file.path(
    survival_dir,
    paste0(
      survival_gene,
      "_group_counts.csv"
    )
  ),
  row.names = FALSE
)

survival_formula <- as.formula(
  paste0(
    "Surv(",
    survival_time_col,
    ", ",
    survival_status_col,
    ") ~ `",
    survival_gene,
    "`"
  )
)

survival_plot <- ggsurvplot(
  fit = survival_fit,
  data = survival_plot_data,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  legend.title = survival_gene,
  legend.labs = c(
    "Low",
    "High"
  ),
  palette = c(
    "#2C7BB6",
    "#D7191C"
  ),
  xlab = paste(
    "Overall survival",
    survival_output_unit
  ),
  ylab = "Overall survival probability",
  ggtheme = theme_classic(
    base_size = 13
  )
)

pdf(
  file = file.path(
    survival_dir,
    paste0(
      survival_gene,
      "_survival_plot.pdf"
    )
  ),
  width = 7,
  height = 7,
  useDingbats = FALSE
)

print(survival_plot)

dev.off()

