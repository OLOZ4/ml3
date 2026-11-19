# =============================================================================
#  UŽDUOTIS 2: mažinti matmenis — PCA ir UMAP (2 bei 6 dimensijos)
#
#  Iš apdoroto duomenų rinkinio sukuriami 4 nauji rinkiniai:
#    a05_pca_2.csv, a05_pca_6.csv, a05_umap_2.csv, a05_umap_6.csv
#  (išvestis į datasets/). Požymiai matmenų mažinimui: visi skaitiniai
#  stulpeliai; class tik pridedama prie rezultato, niekada nedalyvauja.
#
#  Pagrindinė funkcija: reduce_dimensions(data, output_dir)
#  Posistemės: run_pca(), run_umap(), save_reduced()
# =============================================================================

# ---- šio failo nustatymai (keisk čia) ---------------------------------------
settings_reduce <- list(
  verbose = TRUE,               # TRUE = spausdinti šios užduoties žinutes

  # PCA parametrai
  pca = list(
    center = TRUE,              # atimti požymių vidurkius
    scale  = FALSE              # duomenys jau z-skaliuoti, papildoma nereikia
  ),

  # UMAP parametrai (visi umap::umap.config laukai)
  umap = list(
    n_neighbors = 15,           # kaimynų skaičius (vietinė struktūra)
    min_dist    = 0.1,          # minimalus atstumas tarp taškų
    metric      = "euclidean", # atstumo mata
    seed        = 42            # atkartojamumas
  ),

  # kokių matmenų rinkinius kurti: pavadinimas -> dimensijų skaičius
  dimensions = list(pca = c(2, 6), umap = c(2, 6))
)

# ---- pagalbinė: spausdinti tik kai verbose = TRUE ----------------------------
say <- function(verbose, ...) {
  if (verbose) cat(...)
}

# ------------------------------------------------------------------------------

# posistemė: PCA su nurodytu dimensijų skaičiumi
run_pca <- function(data, n_components) {
  feature_names <- setdiff(names(data), "class")
  # pricip() grąžina s-dev, iš jų suskaičiuojame paaiškintą dispersiją
  result <- prcomp(data[, feature_names, drop = FALSE],
                   center = settings_reduce$pca$center,
                   scale  = settings_reduce$pca$scale,
                   rank.  = n_components)
  list(
    embedding   = result$x[, 1:n_components, drop = FALSE],
    explained   = result$sdev[1:n_components]^2,          # komponentų dispersijos
    total_var   = sum(result$sdev^2)                       # visos dispersijos suma
  )
}

# posistemė: UMAP su nurodytu dimensijų skaičiumi
run_umap <- function(data, n_components) {
  feature_names <- setdiff(names(data), "class")
  umap_config <- umap::umap.defaults
  umap_config$n_components <- n_components
  umap_config$n_neighbors  <- settings_reduce$umap$n_neighbors
  umap_config$min_dist    <- settings_reduce$umap$min_dist
  umap_config$metric     <- settings_reduce$umap$metric
  umap_config$seed       <- settings_reduce$umap$seed
  umap_config$n_epochs   <- 200          # numatytasis pagal duomenų dydį
  result <- umap::umap(data[, feature_names, drop = FALSE],
                       config = umap_config)
  list(embedding = as.matrix(result$layout))
}

# posistemė: įrašyti sumažintą rinkinį su class stulpeliu
save_reduced <- function(embedding, class, method, n_components,
                         output_dir, verbose) {
  column_names <- paste0(method, 1:ncol(embedding))
  colnames(embedding) <- column_names
  reduced <- data.frame(embedding, class = class)
  file_name <- sprintf("a05_%s_%d.csv", method, n_components)
  file_path <- file.path(output_dir, file_name)
  write.csv(reduced, file_path, row.names = FALSE)
  say(verbose, sprintf("  %s: %d eilutės x %d dimensijos -> %s\n",
                       method, nrow(reduced), n_components, file_name))
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# pagrindinė užduoties 2 funkcija
reduce_dimensions <- function(data, output_dir) {

  verbose <- settings_reduce$verbose

  # patikra: ar požymių užtenka didžiausiai prašomai dimensijų kartai
  feature_names <- setdiff(names(data), "class")
  max_needed <- max(unlist(settings_reduce$dimensions))
  if (length(feature_names) < max_needed)
    stop("Reikia bent ", max_needed, " požymių, turime:", length(feature_names))

  # PCA
  for (n in settings_reduce$dimensions$pca) {
    say(verbose, sprintf("PCA -> %d dimensijos...\n", n))
    pca_result <- run_pca(data, n)
    explained_share <- sum(pca_result$explained) / pca_result$total_var
    say(verbose, sprintf("  paaiškinta dispersija: %.1f%%\n",
                         100 * explained_share))
    save_reduced(pca_result$embedding, data$class, "pca", n,
                 output_dir, verbose)
  }

  # UMAP
  for (n in settings_reduce$dimensions$umap) {
    say(verbose, sprintf("UMAP -> %d dimensijos...\n", n))
    umap_result <- run_umap(data, n)
    save_reduced(umap_result$embedding, data$class, "umap", n,
                 output_dir, verbose)
  }

  invisible(NULL)
}