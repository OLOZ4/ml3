config <- list(
  input_path        = "a05/A05.csv",     # žali duomenų failas
  output_dir        = "datasets",        # išvesties aplankas
  output_file       = "a05.csv",         # išvesties failas
  outlier_column    = "price_per_carat", # pagal ką šalinamos išskirtys
  iqr_multiplier    = 1.5,               # IQR daugiklis (Q3 + 1.5 x IQR)
  correlation_limit = 0.95               # požymių porų |r| riba
)

source("src/01_pre_processing.r")
source("src/02_reduce_dimensions.r")
source("src/03_find_optimal_clusters.r")
source("src/04_cluster_datasets.r")
source("src/05_stability_experiment.r")


cat("=== 1.1: duomenų nuskaitymas ===\n")
diamonds <- read_raw_data(config$input_path)

cat("\n=== 1.2: netinkamų reikšmių užpildymas ===\n")
diamonds <- fill_invalid_values(diamonds)

cat("\n=== 1.3: išskirčių šalinimas ===\n")
diamonds <- remove_outliers(diamonds, column = config$outlier_column, iqr_multiplier = config$iqr_multiplier)

cat("\n=== 1.4: koreliuotų požymių šalinimas ===\n")
diamonds <- remove_correlated_features(diamonds, threshold = config$correlation_limit)

cat("\n=== 1.5: Normavimas ===\n")
diamonds <- normalize_z_score(diamonds)

# ------------------------- IŠVESTIS ------------------------------------------
dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
output_path <- file.path(config$output_dir, config$output_file)
write.csv(diamonds, output_path, row.names = FALSE)

cat("\n=== 2: dimenciju mažinimas (PCA / UMAP) ===\n")
reduce_dimensions(diamonds, config$output_dir)

cat("\n=== 3: optimalaus klasterių skaičiaus paieška ===\n")
cluster_estimates <- find_optimal_clusters(diamonds)

cat("\n=== 4: klasterizavimas (k-means / hierarchinis / DBSCAN) ===\n")
clustering_results <- cluster_all_datasets()

cat("\n=== 5: stabilumo eksperimentas ===\n")
stability_results <- run_stability_experiment()

cat("\n============================================================\n")
cat("Baigta!\n")
cat("  Eilutės:  ", nrow(diamonds), "\n")
cat("  Stulpeliai:", paste(names(diamonds), collapse = ", "), "\n")
cat("  Išvestis:  ", output_path, "\n")
cat("============================================================\n")