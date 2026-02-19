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
  expr_matrix <- as.matrix(x[, 7:ncol(x)])
  sample_column_names <- colnames(expr_matrix)
  sample_names <- sub(col_to_sample_regex, "\\1", sample_column_names)
  colnames(expr_matrix) <- sample_names
  rownames(expr_matrix) <- x$Geneid

  metadata <- purrr::map_dfr(sample_names, metadata_mapper)
  rownames(metadata) <- metadata$Sample

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


