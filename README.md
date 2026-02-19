### R Bioinformatics Utils

##### process_counts_table
- Uses fread to load the counts txt file produced by featureCounts into a table.
- Accepts a mapper function as an argument to convert sample names into metadata (Useful when sample names contain batch, treatment, etc... information).
- Returns expression matrix and metadata within a SummarizedExperiemnt object to be used downstream.
