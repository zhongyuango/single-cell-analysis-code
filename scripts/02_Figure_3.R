library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(monocle)
library(tradeSeq)
library(BiocParallel)
library(Matrix)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggridges)
library(patchwork)
library(RColorBrewer)
library(viridis)

set.seed(1234)

project_dir <- "."
input_dir   <- file.path(project_dir, "data")
result_dir  <- file.path(project_dir, "results", "trajectory")
figure_dir  <- file.path(result_dir, "figures")
table_dir   <- file.path(result_dir, "tables")
object_dir  <- file.path(result_dir, "objects")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

seurat_file <- file.path(input_dir, "T_cell_processed.rds")
seurat_obj <- readRDS(seurat_file)

celltype_col <- "manu_new_celltype"
sample_col   <- "orig.ident"
group_col    <- "Type"

if (!group_col %in% colnames(seurat_obj@meta.data)) {
  sample_parts <- str_split_fixed(
    seurat_obj@meta.data[[sample_col]],
    "_",
    2
  )
  seurat_obj@meta.data[[group_col]] <- sample_parts[, 2]
}

seurat_obj$celltype <- seurat_obj@meta.data[[celltype_col]]

# Slingshot
slingshot_reduction <- "umap"
start_cluster <- "Naive_CD8"
approx_points <- 150
selected_lineage <- 4

# Monocle2
monocle_max_cells <- 40000
ordering_qval <- 0.05
pseudotime_qval <- 0.01
num_gene_clusters <- 4
root_state <- NULL

save_pdf <- function(plot, filename, width = 7, height = 6) {
  ggsave(
    filename,
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

add_slingshot_pseudotime <- function(seurat_object, sling_object) {
  pt <- as.data.frame(slingPseudotime(sling_object))
  pt <- pt[colnames(seurat_object), , drop = FALSE]

  for (lineage_name in colnames(pt)) {
    seurat_object[[lineage_name]] <- pt[[lineage_name]]
  }

  return(seurat_object)
}

plot_one_lineage <- function(sling_object, lineage_index, output_file) {
  pt <- slingPseudotime(sling_object)[, lineage_index]
  umap <- reducedDims(sling_object)$UMAP

  point_colors <- rep("grey85", length(pt))
  valid <- !is.na(pt)

  point_colors[valid] <- viridis(100)[
    cut(pt[valid], breaks = 100, include.lowest = TRUE)
  ]

  pdf(output_file, width = 6, height = 5, useDingbats = FALSE)

  plot(
    umap,
    col = point_colors,
    pch = 16,
    cex = 0.5,
    asp = 1,
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    main = paste0("Lineage ", lineage_index)
  )

  lines(
    SlingshotDataSet(sling_object),
    linInd = lineage_index,
    lwd = 2,
    col = "black"
  )

  dev.off()
}

plot_density <- function(seurat_object,
                         lineage_column,
                         group_column,
                         selected_groups = NULL,
                         colors = NULL) {
  df <- data.frame(
    pseudotime = seurat_object@meta.data[[lineage_column]],
    group = seurat_object@meta.data[[group_column]]
  ) |>
    filter(!is.na(pseudotime), !is.na(group))

  if (!is.null(selected_groups)) {
    df <- df |>
      filter(group %in% selected_groups)
  }

  p <- ggplot(
    df,
    aes(
      x = pseudotime,
      fill = group,
      color = group
    )
  ) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    theme_classic() +
    labs(x = "Pseudotime", y = "Density")

  if (!is.null(colors)) {
    p <- p +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors)
  }

  return(p)
}

plot_ridge <- function(seurat_object,
                       lineage_column,
                       celltype_column,
                       facet_column = NULL,
                       min_cells = 2) {
  df <- data.frame(
    pseudotime = seurat_object@meta.data[[lineage_column]],
    celltype = seurat_object@meta.data[[celltype_column]]
  )

  if (!is.null(facet_column)) {
    df$facet_group <- seurat_object@meta.data[[facet_column]]
  }

  df <- df |>
    filter(!is.na(pseudotime), !is.na(celltype))

  if (!is.null(facet_column)) {
    df <- df |>
      group_by(facet_group, celltype) |>
      filter(n() >= min_cells) |>
      ungroup()
  } else {
    df <- df |>
      group_by(celltype) |>
      filter(n() >= min_cells) |>
      ungroup()
  }

  p <- ggplot(
    df,
    aes(
      x = pseudotime,
      y = celltype,
      fill = celltype
    )
  ) +
    geom_density_ridges(
      alpha = 0.7,
      scale = 1.1,
      rel_min_height = 0.01
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    ) +
    labs(x = "Pseudotime", y = NULL)

  if (!is.null(facet_column)) {
    p <- p + facet_wrap(~facet_group, ncol = 1)
  }

  return(p)
}

plot_gene_trends <- function(seurat_object,
                             genes,
                             lineage_column,
                             group_column = NULL,
                             group_colors = NULL) {
  expr <- GetAssayData(
    seurat_object,
    assay = "RNA",
    slot = "data"
  )

  valid_genes <- intersect(genes, rownames(expr))
  plots <- list()

  for (gene in valid_genes) {
    df <- data.frame(
      pseudotime = seurat_object@meta.data[[lineage_column]],
      expression = as.numeric(expr[gene, ])
    )

    if (!is.null(group_column)) {
      df$group <- seurat_object@meta.data[[group_column]]
    }

    df <- df |>
      filter(!is.na(pseudotime))

    if (is.null(group_column)) {
      p <- ggplot(df, aes(pseudotime, expression)) +
        geom_smooth(se = FALSE, linewidth = 1)
    } else {
      p <- ggplot(
        df,
        aes(
          pseudotime,
          expression,
          color = group
        )
      ) +
        geom_smooth(se = FALSE, linewidth = 1)

      if (!is.null(group_colors)) {
        p <- p +
          scale_color_manual(values = group_colors)
      }
    }

    plots[[gene]] <- p +
      theme_classic() +
      ggtitle(gene) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      )
  }

  return(plots)
}

sce <- as.SingleCellExperiment(
  seurat_obj,
  assay = "RNA"
)

reducedDims(sce)$UMAP <- Embeddings(
  seurat_obj,
  reduction = slingshot_reduction
)

colData(sce)$celltype <- seurat_obj$celltype
colData(sce)$orig.ident <- seurat_obj@meta.data[[sample_col]]
colData(sce)$group <- seurat_obj@meta.data[[group_col]]

sce_slingshot <- slingshot(
  sce,
  reducedDim = "UMAP",
  clusterLabels = "celltype",
  start.clus = start_cluster,
  approx_points = approx_points
)

lineage_list <- slingLineages(sce_slingshot)
n_lineages <- ncol(slingPseudotime(sce_slingshot))

writeLines(
  capture.output(lineage_list),
  file.path(table_dir, "Slingshot_lineages.txt")
)


celltypes <- sort(unique(sce_slingshot$celltype))

celltype_colors <- setNames(
  colorRampPalette(
    brewer.pal(8, "Set2")
  )(length(celltypes)),
  celltypes
)

pdf(
  file.path(
    figure_dir,
    "Slingshot_trajectory_by_celltype.pdf"
  ),
  width = 7,
  height = 6,
  useDingbats = FALSE
)

plot(
  reducedDims(sce_slingshot)$UMAP,
  col = celltype_colors[sce_slingshot$celltype],
  pch = 16,
  cex = 0.6,
  asp = 1,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

lines(
  SlingshotDataSet(sce_slingshot),
  lwd = 2,
  col = "black"
)

legend(
  "topright",
  legend = names(celltype_colors),
  col = celltype_colors,
  pch = 16,
  bty = "n",
  cex = 0.6
)

dev.off()

for (lineage_index in seq_len(n_lineages)) {
  plot_one_lineage(
    sce_slingshot,
    lineage_index,
    file.path(
      figure_dir,
      paste0(
        "Slingshot_Lineage",
        lineage_index,
        "_pseudotime.pdf"
      )
    )
  )
}

seurat_obj <- add_slingshot_pseudotime(
  seurat_obj,
  sce_slingshot
)

saveRDS(
  seurat_obj,
  file.path(
    object_dir,
    "T_cell_with_Slingshot_pseudotime.rds"
  )
)

saveRDS(
  sce_slingshot,
  file.path(
    object_dir,
    "T_cell_Slingshot_SCE.rds"
  )
)

selected_lineage_name <- colnames(
  slingPseudotime(sce_slingshot)
)[selected_lineage]


group_colors <- c(
  "Naive" = "#8DD3C7",
  "Immunotherapy" = "#DC050C"
)

p_ridge <- plot_ridge(
  seurat_obj,
  selected_lineage_name,
  "celltype",
  facet_column = group_col
)

save_pdf(
  p_ridge,
  file.path(
    figure_dir,
    paste0(
      selected_lineage_name,
      "_ridge_by_celltype_and_group.pdf"
    )
  ),
  width = 10,
  height = 12
)

p_group <- plot_density(
  seurat_obj,
  selected_lineage_name,
  group_col,
  colors = group_colors
)

save_pdf(
  p_group,
  file.path(
    figure_dir,
    paste0(
      selected_lineage_name,
      "_density_by_group.pdf"
    )
  )
)

selected_celltypes <- c(
  "CTL_CD8",
  "Tex(CXCL13+)_CD8",
  "Tfh",
  "Treg"
)

selected_colors <- setNames(
  colorRampPalette(
    brewer.pal(8, "Set2")
  )(length(selected_celltypes)),
  selected_celltypes
)

p_celltype <- plot_density(
  seurat_obj,
  selected_lineage_name,
  "celltype",
  selected_groups = selected_celltypes,
  colors = selected_colors
)

combined_density <- p_celltype / p_group

save_pdf(
  combined_density,
  file.path(
    figure_dir,
    paste0(
      selected_lineage_name,
      "_selected_celltypes_and_group.pdf"
    )
  ),
  width = 7,
  height = 7
)

genes_plot <- c(
  "CXCL13",
  "PDCD1",
  "BCL6",
  "GZMB",
  "IFNG",
  "TIGIT",
  "LAG3",
  "CCR7",
  "CCL5"
)

gene_plots <- plot_gene_trends(
  seurat_obj,
  genes_plot,
  selected_lineage_name,
  group_col,
  group_colors
)

if (length(gene_plots) > 0) {
  p_genes <- wrap_plots(gene_plots, ncol = 3)

  save_pdf(
    p_genes,
    file.path(
      figure_dir,
      paste0(
        selected_lineage_name,
        "_selected_genes_by_group.pdf"
      )
    ),
    width = 10,
    height = 8
  )
}

counts_matrix <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  slot = "counts"
)

pseudotime_matrix <- slingPseudotime(
  sce_slingshot,
  na = FALSE
)

cell_weights <- slingCurveWeights(
  sce_slingshot
)

common_cells <- Reduce(
  intersect,
  list(
    colnames(counts_matrix),
    rownames(pseudotime_matrix),
    rownames(cell_weights)
  )
)

counts_matrix <- counts_matrix[
  ,
  common_cells,
  drop = FALSE
]

pseudotime_matrix <- pseudotime_matrix[
  common_cells,
  ,
  drop = FALSE
]

cell_weights <- cell_weights[
  common_cells,
  ,
  drop = FALSE
]

hvg <- VariableFeatures(seurat_obj)

if (length(hvg) == 0) {
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    nfeatures = 2000
  )
  hvg <- VariableFeatures(seurat_obj)
}

trade_seq_genes <- intersect(
  head(hvg, 2000),
  rownames(counts_matrix)
)

counts_trade_seq <- counts_matrix[
  trade_seq_genes,
  ,
  drop = FALSE
]

if (.Platform$OS.type == "windows") {
  bpparam <- SnowParam(workers = 4, type = "SOCK")
} else {
  bpparam <- MulticoreParam(workers = 4)
}

trade_seq_fit <- fitGAM(
  counts = counts_trade_seq,
  pseudotime = pseudotime_matrix,
  cellWeights = cell_weights,
  nknots = 5,
  verbose = TRUE,
  BPPARAM = bpparam
)

association_result <- associationTest(
  trade_seq_fit,
  lineages = TRUE,
  l2fc = log2(2)
)

association_result$gene <- rownames(
  association_result
)

association_result <- association_result |>
  arrange(pvalue)

write.csv(
  association_result,
  file.path(
    table_dir,
    "tradeSeq_associationTest.csv"
  ),
  row.names = FALSE
)

saveRDS(
  trade_seq_fit,
  file.path(
    object_dir,
    "tradeSeq_fitGAM.rds"
  )
)


monocle_seurat <- seurat_obj

if (!is.null(monocle_max_cells) &&
    ncol(monocle_seurat) > monocle_max_cells) {
  selected_cells <- sample(
    colnames(monocle_seurat),
    monocle_max_cells
  )

  monocle_seurat <- subset(
    monocle_seurat,
    cells = selected_cells
  )
}

counts_monocle <- GetAssayData(
  monocle_seurat,
  assay = "RNA",
  slot = "counts"
)

counts_monocle <- as(
  counts_monocle,
  "sparseMatrix"
)

pheno_data <- monocle_seurat@meta.data

pd <- new(
  "AnnotatedDataFrame",
  data = pheno_data
)

feature_data <- data.frame(
  gene_short_name = rownames(counts_monocle),
  row.names = rownames(counts_monocle)
)

fd <- new(
  "AnnotatedDataFrame",
  data = feature_data
)

monocle_cds <- newCellDataSet(
  counts_monocle,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

nonzero_genes <- rowSums(
  exprs(monocle_cds)
) > 0

monocle_cds <- monocle_cds[
  nonzero_genes,
  ,
  drop = FALSE
]

gene_total_counts <- rowSums(
  exprs(monocle_cds)
)

monocle_cds <- monocle_cds[
  gene_total_counts >= 10,
  ,
  drop = FALSE
]

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

ordering_formula <- paste0("~", celltype_col)

ordering_test <- differentialGeneTest(
  monocle_cds,
  fullModelFormulaStr = ordering_formula,
  cores = 4
)

ordering_test <- ordering_test |>
  arrange(qval)

write.csv(
  ordering_test,
  file.path(
    table_dir,
    "Monocle2_ordering_gene_test.csv"
  )
)

ordering_genes <- rownames(ordering_test)[
  !is.na(ordering_test$qval) &
    ordering_test$qval < ordering_qval
]

if (length(ordering_genes) < 10) {
  ordering_genes <- head(
    rownames(ordering_test),
    2000
  )
}

monocle_cds <- setOrderingFilter(
  monocle_cds,
  ordering_genes
)

p_ordering <- plot_ordering_genes(
  monocle_cds
)

save_pdf(
  p_ordering,
  file.path(
    figure_dir,
    "Monocle2_ordering_genes.pdf"
  )
)

monocle_cds <- reduceDimension(
  monocle_cds,
  max_components = 2,
  reduction_method = "DDRTree"
)

if (is.null(root_state)) {
  monocle_cds <- orderCells(monocle_cds)
} else {
  monocle_cds <- orderCells(
    monocle_cds,
    root_state = root_state
  )
}

saveRDS(
  monocle_cds,
  file.path(
    object_dir,
    "Monocle2_cds_ordered.rds"
  )
)

write.csv(
  pData(monocle_cds),
  file.path(
    table_dir,
    "Monocle2_cell_metadata_with_pseudotime.csv"
  )
)

p_celltype <- plot_cell_trajectory(
  monocle_cds,
  color_by = celltype_col,
  cell_size = 1
)

save_pdf(
  p_celltype,
  file.path(
    figure_dir,
    "Monocle2_trajectory_by_celltype.pdf"
  ),
  width = 8,
  height = 6
)

p_pseudotime <- plot_cell_trajectory(
  monocle_cds,
  color_by = "Pseudotime",
  cell_size = 1
)

save_pdf(
  p_pseudotime,
  file.path(
    figure_dir,
    "Monocle2_trajectory_by_pseudotime.pdf"
  )
)

p_state <- plot_cell_trajectory(
  monocle_cds,
  color_by = "State",
  cell_size = 1
)

save_pdf(
  p_state,
  file.path(
    figure_dir,
    "Monocle2_trajectory_by_state.pdf"
  )
)

if (group_col %in% colnames(pData(monocle_cds))) {
  p_group_monocle <- plot_cell_trajectory(
    monocle_cds,
    color_by = group_col,
    cell_size = 1
  )

  save_pdf(
    p_group_monocle,
    file.path(
      figure_dir,
      "Monocle2_trajectory_by_group.pdf"
    )
  )
}

state_celltype_table <- table(
  pData(monocle_cds)$State,
  pData(monocle_cds)[[celltype_col]]
)

write.csv(
  state_celltype_table,
  file.path(
    table_dir,
    "Monocle2_State_by_celltype.csv"
  )
)

pseudotime_gene_test <- differentialGeneTest(
  monocle_cds,
  fullModelFormulaStr = "~sm.ns(Pseudotime)",
  cores = 4
)

pseudotime_gene_test <- pseudotime_gene_test |>
  arrange(qval)

write.csv(
  pseudotime_gene_test,
  file.path(
    table_dir,
    "Monocle2_pseudotime_gene_test.csv"
  )
)

pseudotime_genes <- rownames(
  pseudotime_gene_test
)[
  !is.na(pseudotime_gene_test$qval) &
    pseudotime_gene_test$qval <
      pseudotime_qval
]

write.table(
  pseudotime_genes,
  file.path(
    table_dir,
    "Monocle2_significant_pseudotime_genes.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

heatmap_genes <- head(
  pseudotime_genes,
  1000
)

if (length(heatmap_genes) >= 2) {
  heatmap_result <- plot_pseudotime_heatmap(
    monocle_cds[heatmap_genes, ],
    num_clusters = num_gene_clusters,
    show_rownames = FALSE,
    return_heatmap = TRUE,
    cluster_rows = TRUE,
    hmcols = colorRampPalette(
      rev(
        brewer.pal(
          n = 11,
          name = "RdYlGn"
        )
      )
    )(50)
  )

  pdf(
    file.path(
      figure_dir,
      "Monocle2_pseudotime_gene_heatmap.pdf"
    ),
    width = 7,
    height = 10,
    useDingbats = FALSE
  )

  print(heatmap_result)
  dev.off()

  gene_clusters <- cutree(
    heatmap_result$tree_row,
    k = num_gene_clusters
  )

  gene_cluster_table <- data.frame(
    gene = names(gene_clusters),
    gene_cluster = as.integer(gene_clusters)
  )

  write.csv(
    gene_cluster_table,
    file.path(
      table_dir,
      "Monocle2_pseudotime_gene_clusters.csv"
    ),
    row.names = FALSE
  )
}

monocle_plot_data <- pData(monocle_cds)

p_monocle_ridge <- ggplot(
  monocle_plot_data,
  aes(
    x = Pseudotime,
    y = .data[[celltype_col]],
    fill = .data[[celltype_col]]
  )
) +
  geom_density_ridges(
    scale = 1,
    rel_min_height = 0.01
  ) +
  scale_y_discrete(position = "right") +
  scale_x_continuous(position = "top") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )

save_pdf(
  p_monocle_ridge,
  file.path(
    figure_dir,
    "Monocle2_pseudotime_ridge_by_celltype.pdf"
  ),
  width = 7,
  height = 8
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    result_dir,
    "trajectory_sessionInfo.txt"
  )
)