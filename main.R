# Generate test data
gwas_data <- generate_gwas_data_for_plot(n_snps = 10000, seed = 123)

# Check the structure
str(gwas_data)

# Test with the Manhattan plot function
plot_manhattan(
  gwas_data = gwas_data,
  chr_col = "CHR",
  bp_col = "BP",
  p_col = "P",
  gene_col = "GENE",
  n_cases = 5000,
  n_controls = 5000,
  plot_title_prefix = "Test GWAS",
  file_name_prefix = "test_gwas"
)

colnames(gwasdataseta)<-c("CHR","BP","P","GENE")



source("R/new_cluster_snps.R")
head(gwasdataseta)
#renv::install("igraph")
library(igraph)
clusterSNPs<-cluster_snps(data = gwasdataseta, rsid_col ="rsid" ,chr_col = "CHR", pos_col = "BP",pvalue_col = "P", pvalue_threshold = 5E-5,distance_threshold = 300000 )

head(clusterSNPs)
create_inverted_manhattan_pair(clusterSNPs, clusterSNPs,
                                           chr_col = "CHR",
                                           bp_col = "BP",
                                           p_col = "P",
                                           gene_col = "GENE",
                                           group_col = "cluster",
                                           loci_to_label = NULL,
                                           genes_to_label = NULL,
                                           custom_gene_colors = NULL,
                                           n_cases1=1000,
                                           n_controls1=1000,
                                           n_cases2=1000,
                                           n_controls2=1000,
                                           plot_title1 = "Dataset 1",
                                           plot_title2 = "Dataset 2",
                                           lambda1 = NULL, lambda2 = NULL,
                                           total_snps_in_study1 = NULL, total_snps_in_study2 = NULL,
                                           add_date_to_title = FALSE,
                                           plot_pval_threshold = 1,
                                           sig_lines = c("red" = 5e-8, "blue" = 1e-5),
                                           label_threshold_colors = c("red" = 5e-8, "orange" = 1e-6, "darkblue" = 1e-5),
                                           y_axis_squish_threshold = 30,
                                           gene_label_size = 3.5,
                                           output_folder = "Inverted_Manhattan_Plots",
                                           file_name_prefix = "inverted_manhattan",
                                           font_family = "Arial Unicode MS")


## Read GWAS dataset from sampledata/gwasdataseta.txt and add to package data
gwasdataseta <- read.table("sampledata/gwasdataseta.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# the chr should be character
gwasdataseta$chr <- as.character(gwasdataseta$chr)
# add random gene names for each 10000 for each SNP
gwasdataseta$gene <- paste0("Gene", seq_len(nrow(gwasdataseta)))
# add random rsid
gwasdataseta$rsid <- paste0("rs", seq_len(nrow(gwasdataseta)))
# Save as package data using usethis (run this script interactively)
if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(gwasdataseta, overwrite = TRUE)
}

# Check structure
str(gwasdataseta)

# Check if 'gwasdataseta' is available in the package data
if ("gwasdataseta" %in% data(package = "manhattwin")$results[, "Item"]) {
  message("gwasdataseta is available in the package data.")
} else {
  warning("gwasdataseta is NOT available in the package data.")
}






unique(clusterSNPs$cluster)

clusterSNPs

