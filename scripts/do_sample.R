# *************************************************************************
# Function to run the sampling procedure
# *************************************************************************

# Parameters:
# sp_code   - species codes column
# min_area  - minimum sample area value
# plot_area - plot area column
# sum_abund - total plot abundance
# repl      - number of replications; defaults to 1000 if not specified

# Returns:
# a 3 element list with sampled abundance mean and SD for each species


do_sample <- function(sp_code, min_area, plot_area, sum_abund, repl = 1000) {
  # Compute sample size
  sample_size <- round(min_area / plot_area[1] * sum_abund[1])
  
  # Sample n times the species names. Will get a matrix of species names.
  # Note: function `replicate()` will simplify its output to matrix or vector.
  # It simplifies to vector if there is one row matrix!
  # This happens for one-species plots.
  # Therefore, the solution with `simplify = FALSE` fallowed by a `do.call`
  # So, operating on lists is safer - https://stackoverflow.com/a/14490315/5193830
  set.seed(666)
  tbl_sp <- replicate(
    n = repl,
    expr = sample(x = sp_code,
                  size = sample_size,
                  replace = TRUE),
    simplify = FALSE
  )
  tbl_sp <- do.call(cbind, tbl_sp)
  
  # Get counts by species from the matrix of species names
  sp_unq <- sort(unique(sp_code))
  tbl_counts <- apply(X = tbl_sp, 
                      MARGIN = 2, 
                      FUN = function(col) table(factor(col, levels = sp_unq)))
  # table(factor(col, levels...)) inspired from
  # https://stat.ethz.ch/pipermail/r-help/2010-June/242381.html
  # Note: function `apply` suffers of the same simplify issue, fix with an if;
  # if vector (does not have dimensions) then convert to one row matrix.
  tbl_counts <- if (is.null(dim(tbl_counts))) matrix(tbl_counts, nrow = 1) else tbl_counts
  
  # Compute sampled mean abundance and SD.
  # Wrap the results with list() so that they pass as columns in the data.table
  return(
    list(sp_code = sp_unq, 
         sampled_abund_mean = matrixStats::rowMeans2(tbl_counts), 
         sampled_abund_sd = matrixStats::rowSds(tbl_counts))
  )
}
