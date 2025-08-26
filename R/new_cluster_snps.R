cluster_snps <- function(data, chr_col = "chr", pos_col = "pos", rsid_col = "rsid",
                         pvalue_col = "pvalue", gene_col = "gene",
                         pvalue_threshold = 0.05, distance_threshold = 1000000) {

  # Load required packages
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required. Please install it with: install.packages('igraph')")
  }

  # Create a copy of the original data to preserve all rows
  result_data <- data

  # Initialize cluster and top columns
  result_data$cluster <- NA_character_
  result_data$top <- NA_character_

  # Check if rsid column exists, if not create a temporary one
  if (!rsid_col %in% names(data)) {
    # Create temporary rsid using chr:pos
    temp_rsid <- paste0(data[[chr_col]], ":", data[[pos_col]])

    # Check for duplicates in temp_rsid
    if (any(duplicated(temp_rsid))) {
      # For duplicates, add unique suffix
      temp_rsid <- ave(temp_rsid, temp_rsid, FUN = function(x) {
        if (length(x) > 1) {
          paste0(x, "_", seq_along(x))
        } else {
          x
        }
      })
    }

    # Add temporary rsid to data
    data$temp_rsid <- temp_rsid
    rsid_col_internal <- "temp_rsid"
  } else {
    rsid_col_internal <- rsid_col
  }

  # Filter data by p-value threshold
  filtered_data <- data[data[[pvalue_col]] <= pvalue_threshold, ]

  # If no SNPs pass the threshold, return original data with NA columns
  if (nrow(filtered_data) == 0) {
    return(result_data)
  }

  # Get unique chromosomes
  chromosomes <- unique(filtered_data[[chr_col]])

  # Process each chromosome separately
  for (chr in chromosomes) {
    # Get SNPs on this chromosome
    chr_snps <- filtered_data[filtered_data[[chr_col]] == chr, ]

    # Skip if only one SNP (can't form clusters)
    if (nrow(chr_snps) <= 1) {
      # For single SNP, create a cluster with just that SNP
      if (nrow(chr_snps) == 1) {
        cluster_id <- paste0("chr", as.character(chr), "_cluster1")
        idx <- which(result_data[[rsid_col_internal]] == chr_snps[[rsid_col_internal]])
        result_data$cluster[idx] <- cluster_id
        result_data$top[idx] <- "top"
      }
      next
    }

    # Create a distance matrix
    positions <- chr_snps[[pos_col]]
    dist_matrix <- outer(positions, positions, function(x, y) abs(x - y))

    # Create adjacency matrix (1 if distance <= threshold, 0 otherwise)
    adj_matrix <- ifelse(dist_matrix <= distance_threshold, 1, 0)
    diag(adj_matrix) <- 0  # Remove self-loops

    # Create graph from adjacency matrix
    g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

    # Find connected components (clusters)
    clusters <- igraph::components(g)

    # For each cluster, identify the top SNP (smallest p-value)
    for (i in 1:clusters$no) {
      # Get indices of SNPs in this cluster
      cluster_indices <- which(clusters$membership == i)

      # Get the SNP with the smallest p-value in this cluster
      top_idx <- cluster_indices[which.min(chr_snps[[pvalue_col]][cluster_indices])]

      # Create cluster ID
      cluster_id <- paste0("chr", as.character(chr), "_cluster", i)

      # Update cluster assignments for all SNPs in this cluster
      for (idx in cluster_indices) {
        rsid <- chr_snps[[rsid_col_internal]][idx]
        result_idx <- which(result_data[[rsid_col_internal]] == rsid)
        result_data$cluster[result_idx] <- cluster_id

        # Mark the top SNP
        if (idx == top_idx) {
          result_data$top[result_idx] <- "top"
        }
      }
    }
  }

  # Remove temporary rsid column if it was created
  if ("temp_rsid" %in% names(result_data)) {
    result_data$temp_rsid <- NULL
  }

  return(result_data)
}
