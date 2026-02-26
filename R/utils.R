#' Process txt output of featureCounts from either a data.table or file path
#' @param x Either a data.table or character path
#'
#' @export
process_counts_table <- function(x, ...) {
  UseMethod("process_counts_table")
}

#' @param x Character path to featureCounts output file
#' @param ... Additional arguments passed to data.table method
#' @export
process_counts_table.character <- function(x, ...){
  process_counts_table(data.table::fread(x), ...)
}

#' @param x data.table containing featureCounts output
#' @param col_to_sample_regex Regular expression with capture group to extract sample names from column names
#' @param metadata_mapper Function that takes sample names and returns named list with extracted metadata
#' @param orgdb Optional OrgDb object for gene annotation (e.g., org.Hs.eg.db, org.Mm.eg.db). If NULL, no symbols added.
#' @return SummarizedExperiment object with counts assay and sample metadata
#' @export
process_counts_table.data.table <- function(x,
                                            col_to_sample_regex = ".*/(.*?)_Aligned.*$",
                                            metadata_mapper,
                                            orgdb = NULL){
  expr_matrix <- to_expr_matrix(x, col_to_sample_regex)
  metadata <- perform_metadata_mapping(colnames(expr_matrix), metadata_mapper)

  gene_info <- x[, 1:6]

  if (!is.null(orgdb)) {
    gene_info$Symbol <- AnnotationDbi::mapIds(
      orgdb,
      keys = x$Geneid,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  }

  rownames(gene_info) <- x$Geneid

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = expr_matrix),
    colData = metadata,
    rowData = gene_info
  )
}

perform_metadata_mapping <- function(sample_names, metadata_mapper){
  metadata <- purrr::map_dfr(sample_names, metadata_mapper)
  rownames(metadata) <- metadata$Sample
  metadata
}

to_expr_matrix <- function(x, col_to_sample_regex){
  expr_matrix <- as.matrix(x[, 7:ncol(x)])
  sample_column_names <- colnames(expr_matrix)
  sample_names <- sub(col_to_sample_regex, "\\1", sample_column_names)
  colnames(expr_matrix) <- sample_names
  rownames(expr_matrix) <- x$Geneid
  expr_matrix
}

#' Label strategies for volcano plots
#'
#' @description
#' Predefined strategies for selecting which genes to label in volcano plots.
#'
#' @format A list of functions:
#' \describe{
#'   \item{candidate_pval}{Label top significant candidates by p-value}
#'   \item{pval_only}{Label top genes by p-value regardless of candidate status}
#'   \item{candidate_lfc}{Label top significant candidates by fold change magnitude}
#' }
#' @export
label_strategies <- list(
  candidate_pval = function(data, n_labels) {
    up_genes <- data |>
      dplyr::filter(candidate, direction == "up") |>
      dplyr::slice_min(padj, n = n_labels, with_ties = FALSE)

    down_genes <- data |>
      dplyr::filter(candidate, direction == "down") |>
      dplyr::slice_min(padj, n = n_labels, with_ties = FALSE)

    dplyr::bind_rows(up_genes, down_genes)
  },

  pval_only = function(data, n_labels) {
    up_genes <- data |>
      dplyr::filter(direction == "up") |>
      dplyr::slice_min(padj, n = n_labels, with_ties = FALSE)

    down_genes <- data |>
      dplyr::filter(direction == "down") |>
      dplyr::slice_min(padj, n = n_labels, with_ties = FALSE)

    dplyr::bind_rows(up_genes, down_genes)
  },

  candidate_lfc = function(data, n_labels) {
    up_genes <- data |>
      dplyr::filter(candidate, direction == "up") |>
      dplyr::slice_max(abs(log2FoldChange), n = n_labels, with_ties = FALSE)

    down_genes <- data |>
      dplyr::filter(candidate, direction == "down") |>
      dplyr::slice_max(abs(log2FoldChange), n = n_labels, with_ties = FALSE)

    dplyr::bind_rows(up_genes, down_genes)
  }
)

#' Create a volcano plot from DESeq2 results
#' @param dds DESeqDataSet object (output from DESeq())
#' @param contrast Optional contrast for results extraction (passed to DESeq2::results())
#' @param lfc_threshold Log2 fold change threshold for significance (default: 1)
#' @param padj_threshold Adjusted p-value threshold for significance (default: 0.05)
#' @param n_labels Number of top genes to label per direction (default: 6)
#' @param label_strategy Function for selecting genes to label. Use label_strategies$candidate_pval (default), label_strategies$pval_only, label_strategies$candidate_lfc, or a custom function
#' @param label_ensembl If TRUE, label genes without symbols using Ensembl IDs. If FALSE, only label genes with symbols (default: FALSE)
#' @param up_color Color for upregulated genes (default: "#FF7777")
#' @param down_color Color for downregulated genes (default: "#7DA8E6")
#' @return A tidyplot/ggplot2 object
#' @export
volcano_plot <- function(dds,
                         contrast = NULL,
                         lfc_threshold = 1,
                         padj_threshold = 0.05,
                         n_labels = 6,
                         label_strategy = label_strategies$candidate_pval,
                         label_ensembl = FALSE,
                         up_color = "#FF7777",
                         down_color = "#7DA8E6") {
  # Extract results
  if (is.null(contrast)) {
    res <- DESeq2::results(dds)
  } else {
    res <- DESeq2::results(dds, contrast = contrast)
  }

  # Convert to data frame and add gene symbols
  df <- as.data.frame(res) |>
    tibble::rownames_to_column("gene_id") |>
    dplyr::mutate(
      neg_log10_padj = -log10(padj),
      direction = dplyr::if_else(log2FoldChange > 0, "up", "down", NA),
      candidate = abs(log2FoldChange) >= lfc_threshold & padj < padj_threshold
    )

  # Add gene symbols if available
  if ("Symbol" %in% colnames(SummarizedExperiment::rowData(dds))) {
    df$Symbol <- SummarizedExperiment::rowData(dds)[df$gene_id, "Symbol"]
  } else {
    df$Symbol <- df$gene_id
  }

  # Create label column based on label_ensembl parameter
  if (label_ensembl) {
    df$Label <- ifelse(is.na(df$Symbol) | df$Symbol == "", df$gene_id, df$Symbol)
  } else {
    df$Label <- df$Symbol
  }

  # Create plot
  df |>
    tidyplots::tidyplot(x = log2FoldChange, y = neg_log10_padj) |>
    tidyplots::add_data_points(
      data = tidyplots::filter_rows(!candidate),
      color = "lightgrey",
      rasterize = TRUE
    ) |>
    tidyplots::add_data_points(
      data = tidyplots::filter_rows(candidate, direction == "up"),
      color = up_color,
      alpha = 0.5
    ) |>
    tidyplots::add_data_points(
      data = tidyplots::filter_rows(candidate, direction == "down"),
      color = down_color,
      alpha = 0.5
    ) |>
    tidyplots::add_reference_lines(
      x = c(-lfc_threshold, lfc_threshold),
      y = -log10(padj_threshold)
    ) |>
    tidyplots::add_data_labels_repel(
      data = function(data) label_strategy(data, n_labels),
      label = Label,
      color = "#000000",
      min.segment.length = 0,
      background = TRUE
    ) |>
    tidyplots::adjust_x_axis_title("$Log[2]~fold~change$") |>
    tidyplots::adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$")
}

#' Create a volcano plot from limma results
#' @param fit MArrayLM object (output from limma::eBayes())
#' @param coef Coefficient or contrast to extract results for (passed to limma::topTable())
#' @param lfc_threshold Log2 fold change threshold for significance (default: 1)
#' @param padj_threshold Adjusted p-value threshold for significance (default: 0.05)
#' @param n_labels Number of top genes to label per direction (default: 6)
#' @param label_strategy Function for selecting genes to label. Use label_strategies$candidate_pval (default), label_strategies$pval_only, label_strategies$candidate_lfc, or a custom function
#' @param up_color Color for upregulated genes (default: "#FF7777")
#' @param down_color Color for downregulated genes (default: "#7DA8E6")
#' @return A tidyplot/ggplot2 object
#' @export
volcano_plot_limma <- function(fit,
                               coef = NULL,
                               lfc_threshold = 1,
                               padj_threshold = 0.05,
                               n_labels = 6,
                               label_strategy = label_strategies$candidate_pval,
                               up_color = "#FF7777",
                               down_color = "#7DA8E6") {
  df <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "P") |>
    tibble::rownames_to_column("gene_id") |>
    dplyr::rename(log2FoldChange = logFC, padj = adj.P.Val) |>
    dplyr::mutate(
      neg_log10_padj = -log10(padj),
      direction = dplyr::if_else(log2FoldChange > 0, "up", "down", NA),
      candidate = abs(log2FoldChange) >= lfc_threshold & padj < padj_threshold,
      Label = gene_id
    )

  df |>
    tidyplots::tidyplot(x = log2FoldChange, y = neg_log10_padj) |>
    tidyplots::add_data_points(
      data = tidyplots::filter_rows(!candidate),
      color = "lightgrey",
      rasterize = TRUE
    ) |>
    tidyplots::add_data_points(
      data = tidyplots::filter_rows(candidate, direction == "up"),
      color = up_color,
      alpha = 0.5
    ) |>
    tidyplots::add_data_points(
      data = tidyplots::filter_rows(candidate, direction == "down"),
      color = down_color,
      alpha = 0.5
    ) |>
    tidyplots::add_reference_lines(
      x = c(-lfc_threshold, lfc_threshold),
      y = -log10(padj_threshold)
    ) |>
    tidyplots::add_data_labels_repel(
      data = function(data) label_strategy(data, n_labels),
      label = Label,
      color = "#000000",
      min.segment.length = 0,
      background = TRUE
    ) |>
    tidyplots::adjust_x_axis_title("$Log[2]~fold~change$") |>
    tidyplots::adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$")
}


#' Plot UMAP with cluster hulls
#'
#' @param umap_df Data frame containing UMAP1, UMAP2, cluster, and value columns
#' @param alpha Alpha shape parameter for hull construction (default: 5)
#' @param buffer_dist Buffer distance to expand hulls (default: 0.3)
#' @return A tidyplot/ggplot2 object
#' @export
plot_umap_with_hulls <- function(umap_df, alpha = 5, buffer_dist = 0.3) {
  order_ahull_edges <- function(edges) {
    find_next <- function(remaining, current, threshold = 1e-10) {
      which(
        (abs(remaining$x1 - current[1]) < threshold & abs(remaining$y1 - current[2]) < threshold) |
          (abs(remaining$x2 - current[1]) < threshold & abs(remaining$y2 - current[2]) < threshold)
      )
    }
    orient_edge <- function(row, current, threshold = 1e-10) {
      if (abs(row$x2 - current[1]) < threshold & abs(row$y2 - current[2]) < threshold)
        row[, c("x1","y1","x2","y2")] <- row[, c("x2","y2","x1","y1")]
      row
    }
    walk_edges <- function(result, remaining, current) {
      idx <- find_next(remaining, current)
      if (length(idx) == 0 || nrow(remaining) == 0) return(result)
      next_edge <- orient_edge(remaining[idx[1], ], current)
      walk_edges(rbind(result, next_edge), remaining[-idx[1], ], c(next_edge$x2, next_edge$y2))
    }
    ordered <- walk_edges(edges[1, ], edges[-1, ], c(edges$x2[1], edges$y2[1]))
    tibble::tibble(UMAP1 = c(ordered$x1, ordered$x1[1]),
                   UMAP2 = c(ordered$y1, ordered$y1[1]))
  }

  expand_hull <- function(hull_df, buffer_dist) {
    hull_df |>
      as.matrix() |>
      list() |>
      sf::st_polygon() |>
      sf::st_buffer(buffer_dist) |>
      sf::st_coordinates() |>
      tibble::as_tibble() |>
      dplyr::select(UMAP1 = X, UMAP2 = Y)
  }

  cluster_hulls <- umap_df |>
    dplyr::group_by(cluster) |>
    dplyr::group_modify(~ {
      ah <- alphahull::ashape(as.matrix(.x[, c("UMAP1", "UMAP2")]), alpha = alpha)
      edges <- as.data.frame(ah$edges[, c("x1", "y1", "x2", "y2")])
      ordered_hull <- order_ahull_edges(edges)
      expand_hull(ordered_hull, buffer_dist = buffer_dist)
    }) |>
    dplyr::ungroup()

  centroids <- umap_df |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(dplyr::across(c(UMAP1, UMAP2), mean))

  umap_df |>
    tidyplots::tidyplot(x = UMAP1, y = UMAP2, color = value) |>
    tidyplots::add(ggnewscale::new_scale_fill()) |>
    tidyplots::add(ggplot2::geom_polygon(
      data = cluster_hulls,
      ggplot2::aes(x = UMAP1, y = UMAP2, group = cluster, fill = cluster),
      alpha = 0.5, color = "black", linewidth = 0.5,
      inherit.aes = FALSE
    )) |>
    tidyplots::add(ggplot2::scale_fill_brewer(palette = "Pastel1", name = "Cluster")) |>
    tidyplots::add_data_points(size = 1.2, alpha = 0.6) |>
    tidyplots::add(ggplot2::geom_label(
      data = centroids,
      ggplot2::aes(x = UMAP1, y = UMAP2, label = cluster),
      inherit.aes = FALSE, size = 3, alpha = 0.95
    ))
}

