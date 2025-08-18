#' Generate Well-Distributed GWAS Data for Plotting
#'
#' Creates a synthetic GWAS dataset specifically designed for generating visually
#' appealing and realistic Manhattan plots. The p-values are generated in tiers
#' to ensure a clear distribution of non-significant, suggestive, and genome-wide
#' significant signals.
#'
#' @param n_snps Total number of SNPs to generate (default: 20000).
#' @param n_hits Number of genome-wide significant "hits" to generate (default: 10).
#' @param n_suggestive Number of suggestive signals to generate (default: 50).
#' @param gws_p The genome-wide significance p-value threshold (default: 5e-8).
#' @param suggestive_p The suggestive significance p-value threshold (default: 1e-5).
#' @param seed Random seed for reproducibility (default: NULL).
#' @param include_x Logical: should chromosome X be included? (default: TRUE).
#'
#' @return A data.frame with CHR, BP, P, and GENE columns, ready for plotting.
#'
#' @examples
#' # Generate a standard dataset with 10 top hits
#' gwas_plot_data <- generate_gwas_data_for_plot(seed = 42)
#'
#' # Generate a denser plot with more suggestive signals
#' gwas_dense_data <- generate_gwas_data_for_plot(n_snps = 50000, n_suggestive = 200, seed = 123)
#'
#' # To create a mirrored plot (as in the example image), you would generate
#' # two datasets and plot one with negative -log10(P) values.
#' dataset1 <- generate_gwas_data_for_plot(seed = 1)
#' dataset2 <- generate_gwas_data_for_plot(seed = 2)
#'
#' @export
generate_gwas_data_for_plot <- function(n_snps = 20000,
                                        n_hits = 10,
                                        n_suggestive = 50,
                                        gws_p = 5e-8,
                                        suggestive_p = 1e-5,
                                        seed = NULL,
                                        include_x = TRUE) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Ensure the number of signals is not more than the total SNPs
  if (n_hits + n_suggestive > n_snps) {
    stop("The sum of n_hits and n_suggestive cannot exceed n_snps.")
  }

  # --- 1. Define Chromosome Structure ---
  chromosomes <- as.character(1:22)
  chr_lengths <- c(
    "1" = 249, "2" = 243, "3" = 198, "4" = 191, "5" = 181, "6" = 171,
    "7" = 159, "8" = 146, "9" = 141, "10" = 135, "11" = 135, "12" = 133,
    "13" = 115, "14" = 107, "15" = 102, "16" = 90, "17" = 81, "18" = 78,
    "19" = 59, "20" = 63, "21" = 48, "22" = 51
  )
  if (include_x) {
    chromosomes <- c(chromosomes, "X")
    chr_lengths["X"] <- 155
  }

  # --- 2. Generate SNP Positions (Vectorized) ---
  chr_probs <- chr_lengths / sum(chr_lengths)
  snp_chr <- sample(chromosomes, size = n_snps, replace = TRUE, prob = chr_probs)

  max_bp_for_snps <- chr_lengths[snp_chr] * 1e6
  snp_bp <- floor(runif(n_snps, min = 1, max = max_bp_for_snps + 1))

  # --- 3. Generate P-values in Tiers ---
  p_values <- rep(1.0, n_snps)
  all_indices <- 1:n_snps

  # Select indices for hits
  hit_indices <- sample(all_indices, size = n_hits)
  remaining_indices <- setdiff(all_indices, hit_indices)

  # Select indices for suggestive SNPs from the remainder
  suggestive_indices <- sample(remaining_indices, size = n_suggestive)

  # The rest are background noise
  background_indices <- setdiff(remaining_indices, suggestive_indices)

  # Generate p-values for each tier
  # For hits, p-values are well below the GWS line
  p_values[hit_indices] <- runif(n_hits, min = 1e-20, max = gws_p * 0.9)

  # For suggestive, p-values are between the suggestive and GWS lines
  p_values[suggestive_indices] <- runif(n_suggestive, min = gws_p, max = suggestive_p)

  # For background, p-values are above the suggestive line
  p_values[background_indices] <- runif(length(background_indices), min = suggestive_p, max = 1)

  # --- 4. Add Gene Labels for Top Hits ---
  genes <- rep(NA_character_, n_snps)
  # Generate unique names like "GENE001", "GENE002", etc.
  hit_gene_names <- paste0("GENE", sprintf(paste0("%0", nchar(n_hits), "d"), 1:n_hits))
  genes[hit_indices] <- hit_gene_names

  # --- 5. Assemble the Final Data Frame ---
  gwas_data <- data.frame(
    CHR = snp_chr,
    BP = snp_bp,
    P = p_values,
    GENE = genes,
    stringsAsFactors = FALSE
  )

  # Sort by chromosome and position, which is standard for GWAS data
  gwas_data$CHR <- factor(gwas_data$CHR, levels = chromosomes)
  gwas_data <- gwas_data[order(gwas_data$CHR, gwas_data$BP), ]

  return(gwas_data)
}
