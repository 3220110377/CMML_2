library(dplyr)
library(tidyr)
library(ggplot2)

plotMarkerPercentageSeurat <- function(seurat_obj,
                                       markers,
                                       groupClusters,
                                       assayName = "RNA") {
  # 1. Identify assays actually present in the Seurat object
  availableAssays <- names(seurat_obj@assays)
  assayName      <- intersect(assayName, availableAssays)
  if (length(assayName) == 0) {
    stop("Specified assay not found in Seurat object. Available assays: ",
         paste(availableAssays, collapse = ", "))
  }
  
  # 2. Build metadata data.frame
  meta <- data.frame(
    cell    = colnames(seurat_obj),
    cluster = as.character(Idents(seurat_obj)),
    stringsAsFactors = FALSE
  )
  # Map cluster IDs to group labels
  meta <- meta %>%
    rowwise() %>%
    mutate(group = {
      hits <- names(groupClusters)[
        vapply(groupClusters, function(ids) as.character(cluster) %in% as.character(ids), logical(1))
      ]
      if (length(hits)) hits[1] else NA_character_
    }) %>%
    ungroup() %>%
    filter(!is.na(group))
  
  # 3. Calculate percentage of cells expressing each marker set for each assay
  df_list <- list()
  for (assay in assayName) {
    # Extract raw counts matrix for this assay
    mat <- GetAssayData(seurat_obj, assay = assay, slot = "counts")
    for (grp in names(markers)) {
      genes <- intersect(markers[[grp]], rownames(mat))
      if (length(genes) == 0) next
      submat    <- mat[genes, meta$cell, drop = FALSE]
      expr_flag <- colSums(submat > 0, na.rm = TRUE) > 0
      
      tmp <- meta %>%
        mutate(expr = expr_flag[cell]) %>%
        group_by(group) %>%
        summarise(
          n_cells = n(),                         # total cells in group
          n_expr  = sum(expr, na.rm = TRUE),     # cells expressing ≥1 marker
          .groups = "drop"
        ) %>%
        mutate(
          pct        = n_expr / n_cells * 100,   # percentage expressing
          marker_set = grp,
          assay      = assay
        )
      
      df_list[[paste(assay, grp, sep = "_")]] <- tmp
    }
  }
  df_all <- bind_rows(df_list)
  
  # 4. Generate barplot with percentage labels
  ggplot(df_all, aes(x = group, y = pct, fill = assay)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", pct)),
              position = position_dodge(width = 0.7),
              vjust = -0,
              size = 3) +
    facet_wrap(~ marker_set, ncol = 1, scales = "free_y") +
    labs(
      x = "Group",
      y = "Percent of cells expressing ≥1 marker",
      fill = "Assay"
    ) +
    theme_bw() +
    theme(
      panel.grid       = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )
}

plotMarkerViolin <- function(seurat_obj,
                             markers,
                             groupClusters,
                             assay = "RNA",
                             slot = "data",
                             ncol = 3) {
  
  # 1. Build metadata with cell and cell type information
  meta <- data.frame(
    Cell    = colnames(seurat_obj),
    Cluster = as.character(Idents(seurat_obj)),
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(CellType = {
      # Map clusters to human-readable cell type labels
      ct <- names(groupClusters)[
        sapply(groupClusters, function(ids) Cluster %in% as.character(ids))
      ]
      if (length(ct)) ct[1] else NA_character_
    }) %>%
    ungroup() %>%
    filter(!is.na(CellType))
  
  # 2. Extract expression matrix for specified assay and slot
  expr_mat <- SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # 3. Filter to markers present in the data
  all_markers <- unique(unlist(markers))
  keep_genes  <- intersect(all_markers, rownames(expr_mat))
  keep_cells  <- meta$Cell
  sub_mat     <- expr_mat[keep_genes, keep_cells, drop = FALSE]
  
  # 4. Convert to long format for ggplot
  df <- as.data.frame(as.matrix(sub_mat)) %>%
    tibble::rownames_to_column(var = "Marker") %>%
    pivot_longer(
      cols      = -Marker,
      names_to  = "Cell",
      values_to = "Expression"
    ) %>%
    inner_join(meta, by = "Cell")
  
  # 5. Create violin plot of expression by cell type
  p <- ggplot(df, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ Marker, scales = "free_y", ncol = ncol) +
    theme_bw() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x     = "Cell Type",
      y     = paste0("Expression (", slot, ")"),
      title = "Marker Expression Across Cell Types"
    )
  return(p)
}
