# Generate test data
gwas_data <- generate_random_gwas_data(n_snps = 10000, seed = 123)

# Check the structure
str(gwas_data)

# Test with the Manhattan plot function
manplot(
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

create_inverted_manhattan_pair(gwas_data, gwas_data,
                                           chr_col = "CHR",
                                           bp_col = "BP",
                                           p_col = "P",
                                           gene_col = "GENE",
                                           group_col = NULL,
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
