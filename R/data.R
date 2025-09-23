#' A SingleCellExperiment object containing both cell type and pseudotime
#' @format A dataset with 10 rows (genes) and 1289 cols (cells)
#' @return The corresponding SingleCellExperiment object
#' @usage data("example_sce")
"example_sce"

#' A sparse matrix with example data
#' @format A dgCMatrix with 100 rows (genes) and 2087 cols (cells)
#' @return The corresponding sparse matrix object 'example_count'
#' @usage data("example_count")
"example_count"

#' A numeric vector of pseudotime values
#' @format A numeric vector with 2087 values
#' @return The corresponding pseudotime of the 'example_count'
#' @usage data("pseudotime")
"pseudotime"