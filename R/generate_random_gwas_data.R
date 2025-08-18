#' Generate Random GWAS Data for Manhattan Plot Testing
#'
#' Creates a synthetic GWAS dataset with chromosome, position, p-value, and optional gene information.
#' Suitable for testing the create_advanced_manhattan function.
#'
#' @param n_snps Number of SNPs to generate (default: 10000)
#' @param chr_col Name for chromosome column (default: "CHR")
#' @param bp_col Name for base pair position column (default: "BP")
#' @param p_col Name for p-value column (default: "P")
#' @param gene_col Name for gene column. Set to NULL to omit (default: "GENE")
#' @param include_x Logical: include chromosome X? (default: TRUE)
#' @param signal_proportion Proportion of SNPs to be significant signals (default: 0.01)
#' @param min_signal_p Minimum p-value for significant signals (default: 1e-10)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return A data.frame with synthetic GWAS data
#'
#' @examples
#' # Generate basic GWAS data
#' gwas_data <- generate_random_gwas_data(n_snps = 5000)
#' 
#' # Generate data without gene information
#' gwas_data_no_genes <- generate_random_gwas_data(gene_col = NULL)
#'
#' @export
generate_random_gwas_data <- function(n_snps = 10000,
                                      chr_col = "CHR",
                                      bp_col = "BP",
                                      p_col = "P",
                                      gene_col = "GENE",
                                      include_x = TRUE,
                                      signal_proportion = 0.01,
                                      min_signal_p = 1e-10,
                                      seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Define chromosomes
  chromosomes <- as.character(1:22)
  if (include_x) {
    chromosomes <- c(chromosomes, "X")
  }
  
  # Create chromosome lengths (approximate in Mb)
  chr_lengths <- c(
    "1" = 248, "2" = 242, "3" = 198, "4" = 190, "5" = 181,
    "6" = 170, "7" = 158, "8" = 146, "9" = 140, "10" = 135,
    "11" = 134, "12" = 132, "13" = 114, "14" = 107, "15" = 101,
    "16" = 90, "17" = 83, "18" = 80, "19" = 58, "20" = 64,
    "21" = 46, "22" = 50, "X" = 154
  )
  
  # Filter chromosomes based on include_x
  chr_lengths <- chr_lengths[names(chr_lengths) %in% chromosomes]
  
  # Distribute SNPs across chromosomes proportional to their length
  chr_probs <- chr_lengths / sum(chr_lengths)
  chr_sample <- sample(names(chr_lengths), size = n_snps, replace = TRUE, prob = chr_probs)
  
  # Generate base pair positions within each chromosome
  bp_positions <- sapply(chr_sample, function(chr) {
    max_bp <- chr_lengths[chr] * 1e6  # Convert Mb to bp
    round(runif(1, min = 1, max = max_bp))
  })
  
  # Determine which SNPs will be signals
  n_signals <- round(n_snps * signal_proportion)
  signal_indices <- sample(seq_len(n_snps), size = n_signals, replace = FALSE)
  
  # Generate p-values
  p_values <- runif(n_snps, min = 0, max = 1)
  
  # Make signal SNPs more significant
  p_values[signal_indices] <- runif(n_signals, min = min_signal_p, max = 5e-8)
  
  # Create the data frame with the correct number of rows
  result <- data.frame(
    matrix(nrow = n_snps, ncol = 0),
    stringsAsFactors = FALSE
  )
  
  # Assign columns using standard data.frame methods
  result[[chr_col]] <- chr_sample
  result[[bp_col]] <- bp_positions
  result[[p_col]] <- p_values
  
  # Add gene information if requested
  if (!is.null(gene_col)) {
    # Sample some genes to label (about 5% of SNPs)
    n_genes <- round(n_snps * 0.05)
    gene_indices <- sample(seq_len(n_snps), size = n_genes, replace = FALSE)
    
    # Generate random gene names
    gene_names <- paste0("GENE", sprintf("%04d", 1:n_genes))
    genes <- rep(NA_character_, n_snps)
    genes[gene_indices] <- sample(gene_names, size = n_genes, replace = TRUE)
    
    # Add some specific gene names for testing custom colors
    if (n_genes >= 5) {
      genes[gene_indices[1:5]] <- c("BRCA1", "APOE", "FTO", "TCF7L2", "PPARG")
    }
    
    result[[gene_col]] <- genes
  }
  
  return(result)
}