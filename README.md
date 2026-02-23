# R Bioinformatics Utils

## Installation
```r
devtools::install("path/to/RBioinfUtils")
```

Requires: `DESeq2`, `SummarizedExperiment`, `data.table`, `purrr`, `dplyr`, `tidyplots`, `alphahull`, `sf`, `ggnewscale`, `ggplot2`, and (optional) organism annotation packages (e.g., `org.Mm.eg.db`, `org.Hs.eg.db`)

## Functions

### `process_counts_table()`

Loads featureCounts output into a `SummarizedExperiment` object for downstream analysis.

- Converts sample names to metadata by using user-provided function that takes a sample name string and returns a named list of metadata (e.g., `Group`, `Batch`, `Replicate`). Useful when sample names encode experimental design information.
- Adds gene symbols using Bioconductor annotation databases
- Returns expression matrix and metadata in a `SummarizedExperiment`

### `volcano_plot()`

Creates volcano plots from DESeq2 results.

- Highlights significant genes by fold change and p-value thresholds
- Labels top genes using customizable strategies
- Supports gene symbols with optional fallback to Ensembl IDs for labeling

**Label strategies:** Control which genes are labeled:
- `label_strategies$candidate_pval` - Top significant genes by p-value (default)
- `label_strategies$candidate_lfc` - Top significant genes by fold change
- `label_strategies$pval_only` - Top genes by p-value (ignores significance)

Custom strategies can be defined as functions that filter and select genes from the data.

### `plot_umap_with_hulls()`

Creates UMAP plots with convex hulls around clusters.

- Draws alpha-shape hulls around cluster point clouds
- Labels cluster centroids
- Colors points by continuous value (e.g., expression, score)

**Required columns in `umap_df`:**
- `UMAP1` - First UMAP dimension
- `UMAP2` - Second UMAP dimension  
- `cluster` - Cluster assignment (factor or character)
- `value` - Continuous variable for point coloring
