# =============================================================================
#  UŽDUOTIS 4: klasterizuoti visus duomenų rinkinius trimis algoritmais
#
#  Kiekvienam rinkiniui (a05, pca_2, pca_6, umap_2, umap_6) taikomi:
#    - k-means                (k = cluster_count)
#    - hierarchinis (Ward)    (k = cluster_count)
#    - DBSCAN                 (eps: nustatytas arba automatinis; klasterių
#                              skaičius gaunasi iš duomenų, k = nesvarstomas)
#
#  Išvestis:
#    - 15 vizualizacijų + bendra 3x5 grafikų matrica  -> visuals/<aplankas>/
#    - visų matų lentelė (vidiniai + išoriniai)       -> CSV
#
#  Matai:
#    Vidiniai:  Silueto plotis, Calinski-Harabasz, Davies-Bouldin
#    Išoriniai (palyginus su class): Jaccard, ARI, NMI
#
#  Pagrindinė funkcija: cluster_all_datasets()
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cluster)
  library(dbscan)
  library(FNN)
  library(gridExtra)
})

if (!exists("say")) say <- function(verbose, ...) if (verbose) cat(...)

# ---- šio failo nustatymai (keisk čia) ----------------------------------------
settings_cluster_all <- list(
  verbose = TRUE,               # TRUE = spausdinti šios užduoties žinutes

  cluster_count = 4,             # KLASTERIŲ SKAIČIUS k-means ir hierarchiniam

  # kuriuos rinkinius klasterizuoti (datasets/ aplanke)
  datasets = c("a05.csv", "a05_pca_2.csv", "a05_pca_6.csv",
               "a05_umap_2.csv", "a05_umap_6.csv"),
  data_dir = "datasets",

  # algoritmų parametrai
  kmeans = list(nstart = 25, seed = 75),
  hierarchical = list(method = "ward.D2"),
  dbscan = list(
    minPts = 5,                 # minimalus taškų skaičius „tankioje“ zonoje
    eps    =0.8              # NULL = automatinis (kNN atstumų kelio taškas)
  ),                            # arba skaičius, pvz. 0.5

  visual_dir = "visuals/clustering_summary",
  metrics_file = "clustering_metrics.csv",
  plot = list(width = 10, height = 8, dpi = 150,         # pavienių grafikų
              matrix_width = 16, matrix_height = 22)     # bendros matricos
)

# ==============================================================================
# POSISTEMĖS: algoritmai
# ==============================================================================

run_kmeans_clustering <- function(features) {
  set.seed(settings_cluster_all$kmeans$seed)
  kmeans(features,
         centers = settings_cluster_all$cluster_count,
         nstart  = settings_cluster_all$kmeans$nstart,
         iter.max = 100)$cluster
}

run_hierarchical_clustering <- function(features) {
  hierarchy <- hclust(dist(features),
                      method = settings_cluster_all$hierarchical$method)
  cutree(hierarchy, k = settings_cluster_all$cluster_count)
}

run_dbscan_clustering <- function(features) {
  min_pts <- settings_cluster_all$dbscan$minPts
  eps     <- settings_cluster_all$dbscan$eps

  # automatinis eps: kNN atstumų (k = minPts) kreivės kelio taškas –
  # taškas, toliausiai nuo tiesės, jungiančios mažiausią ir didžiausią atstumą
  if (is.null(eps)) {
    knn_distances <- sort(get.knn(features, k = min_pts)$nn.dist[, min_pts])
    n <- length(knn_distances)
    slope <- (knn_distances[n] - knn_distances[1]) / (n - 1)
    line_y <- knn_distances[1] + slope * (0:(n - 1))
    eps <- knn_distances[which.max(abs(knn_distances - line_y))]
    say(settings_cluster_all$verbose,
        sprintf("    automatinis eps = %.3f\n", eps))
  }

  result <- dbscan(features, eps = eps, minPts = min_pts)
  as.integer(result$cluster)   # 0 = triukšmas
}

# ==============================================================================
# POSISTEMĖS: vidiniai matai
# ==============================================================================

# Silueto plotis (vidutinis); NA kai klasterių < 2
silhouette_score <- function(features, clusters) {
  if (length(unique(clusters)) < 2) return(NA_real_)
  mean(cluster::silhouette(clusters, dist(features))[, "sil_width"])
}

# Calinski-Harabasz: tarpklasterinės / vidinės dispersijos santykis
calinski_harabasz <- function(features, clusters) {
  cluster_ids <- sort(unique(clusters))   # aggregate grąžina surūšiuotus
  k <- length(cluster_ids)
  n <- nrow(features)
  if (k < 2 || k >= n) return(NA_real_)
  centers <- as.matrix(aggregate(features, by = list(cluster = clusters),
                                  mean)[, -1])
  global_center <- colMeans(features)
  within_ss  <- sum(sapply(seq_len(k), function(i) {
    members <- features[clusters == cluster_ids[i], , drop = FALSE]
    sum((members - centers[i, ])^2)
  }))
  between_ss <- sum(sapply(seq_len(k), function(i) {
    n_i <- sum(clusters == cluster_ids[i])
    n_i * sum((centers[i, ] - global_center)^2)
  }))
  (between_ss / (k - 1)) / (within_ss / (n - k))
}

# Davies-Bouldin: klasterių sklaidos ir centrų atstumo santykis (mažiau = geriau)
davies_bouldin <- function(features, clusters) {
  cluster_ids <- sort(unique(clusters))   # aggregate grąžina surūšiuotus
  k <- length(cluster_ids)
  if (k < 2) return(NA_real_)
  centers <- as.matrix(aggregate(features, by = list(cluster = clusters),
                                  mean)[, -1])
  spread <- sapply(cluster_ids, function(id) {
    members <- features[clusters == id, , drop = FALSE]
    mean(sqrt(rowSums((members - centers[cluster_ids == id, ])^2)))
  })
  pairwise_scores <- sapply(seq_len(k), function(i) {
    max(sapply(setdiff(seq_len(k), i), function(j) {
      (spread[i] + spread[j]) /
        sqrt(sum((centers[i, ] - centers[j, ])^2))
    }))
  })
  mean(pairwise_scores)
}

# ==============================================================================
# POSISTEMĖS: išoriniai matai (klasterizavimas vs žinomos klasės)
# ==============================================================================

# porinės statistikos pagal sutapimo lentelę (klasteriai x klasės)
pair_counts <- function(clusters, labels) {
  contingency <- table(clusters, labels)
  n <- sum(contingency)
  a <- sum(choose(contingency, 2))                      # kartu abejuose
  b <- sum(choose(rowSums(contingency), 2)) - a         # kartu tik klasteriuose
  c <- sum(choose(colSums(contingency), 2)) - a         # kartu tik klasėse
  total <- choose(n, 2)
  d <- total - a - b - c
  list(a = a, b = b, c = c, d = d, total = total)
}

# Jaccard: bendrų elementų dalis (porinis)
jaccard_index <- function(clusters, labels) {
  counts <- pair_counts(clusters, labels)
  counts$a / (counts$a + counts$b + counts$c)
}

# ARI: pakoreguotas Rand indeksas
adjusted_rand_index <- function(clusters, labels) {
  counts <- pair_counts(clusters, labels)
  expected <- (counts$a + counts$b) * (counts$a + counts$c) / counts$total
  max_index <- (counts$a + counts$b + counts$a + counts$c) / 2
  if (max_index == expected) return(NA_real_)
  (counts$a - expected) / (max_index - expected)
}

# NMI: normalizuota tarpusavio informacija
normalized_mutual_info <- function(clusters, labels) {
  contingency <- table(clusters, labels)
  n <- sum(contingency)
  prob_ij <- contingency / n
  prob_i  <- rowSums(contingency) / n
  prob_j  <- colSums(contingency) / n
  mutual_info <- sum(prob_ij *
                     log(ifelse(prob_ij > 0, prob_ij, 1) /
                         (prob_i %*% t(prob_j))))
  entropy_i <- -sum(prob_i * log(ifelse(prob_i > 0, prob_i, 1)))
  entropy_j <- -sum(prob_j * log(ifelse(prob_j > 0, prob_j, 1)))
  if (entropy_i * entropy_j == 0) return(NA_real_)
  mutual_info / sqrt(entropy_i * entropy_j)
}

# ==============================================================================
# POSISTEMĖS: vizualizacija
# ==============================================================================

# vizualinėms plokštumoms: jei rinkinys 2D – jo dvi dimensijos,
# kitaip – pirmos dvi pagrindinės komponentės
viz_coordinates <- function(features) {
  if (ncol(features) == 2) {
    result <- features
  } else {
    result <- prcomp(features)$x[, 1:2]
  }
  data.frame(x = result[, 1], y = result[, 2])
}

# spalvos klasteriams (0 = triukšmas pilkas)
cluster_palette <- function(clusters) {
  unique_ids <- sort(unique(clusters))
  colors <- hcl.colors(max(length(unique_ids), 2), "Dark 3")
  if (0 %in% unique_ids) colors[1] <- "grey60"   # triukšmas pilkas
  setNames(colors[seq_along(unique_ids)], unique_ids)
}

# vienas sklaidos grafikas spalvotas pagal klasterius
plot_cluster_visual <- function(dataset_name, algorithm, viz, clusters,
                               axis_is_pc = TRUE) {
  plot_data <- data.frame(x = viz$x, y = viz$y, cluster = factor(clusters))
  n_clusters <- length(unique(clusters))
  subtitle <- sprintf("%s | klasterių: %d%s", dataset_name, n_clusters,
                      if (any(clusters == 0))
                        sprintf(" (triukšmo: %.1f%%)",
                                100 * mean(clusters == 0)) else "")
  x_label <- if (axis_is_pc) "PK1" else "Dimensija 1"
  y_label <- if (axis_is_pc) "PK2" else "Dimensija 2"
  ggplot(plot_data, aes(x = x, y = y, colour = cluster)) +
    geom_point(size = 1.1, alpha = 0.7) +
    scale_colour_manual(values = cluster_palette(clusters),
                        name = "Klasteris") +
    labs(title = paste0(toupper(algorithm), " – ", dataset_name),
         subtitle = subtitle,
         x = x_label, y = y_label) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "right")
}

# ==============================================================================
# PAGRINDINĖ UŽDUOTIES 4 FUNKCIJA
# ==============================================================================

cluster_all_datasets <- function() {

  verbose <- settings_cluster_all$verbose
  say(verbose, "=== UŽDUOTIS 4: klasterizavimas (k-means / hierarchinis / DBSCAN) ===\n")
  say(verbose, "  klasterių skaičius (k-means ir hierarchiniam):",
      settings_cluster_all$cluster_count, "\n")

  dir.create(settings_cluster_all$visual_dir,
             showWarnings = FALSE, recursive = TRUE)

  algorithms <- list(
    kmeans       = run_kmeans_clustering,
    hierarchical = run_hierarchical_clustering,
    dbscan       = run_dbscan_clustering
  )

  metrics_table <- data.frame()
  all_plots <- list()
  plot_index <- 0

  for (dataset_name in settings_cluster_all$datasets) {

    say(verbose, sprintf("\n  Rinkinys: %s\n", dataset_name))
    dataset <- read.csv(file.path(settings_cluster_all$data_dir, dataset_name),
                        check.names = FALSE)
    feature_names <- setdiff(names(dataset), "class")
    features <- dataset[, feature_names, drop = FALSE]
    labels <- dataset$class
    viz <- viz_coordinates(features)
    ncol <- ncol(features)

    for (algorithm_name in names(algorithms)) {

      say(verbose, sprintf("    %s ...\n", algorithm_name))
      clusters <- algorithms[[algorithm_name]](features)

      # vidiniai matai
      internal <- data.frame(
        silhouette        = silhouette_score(features, clusters),
        calinski_harabasz = calinski_harabasz(features, clusters),
        davies_bouldin    = davies_bouldin(features, clusters)
      )

      # išoriniai matai (palyginus su class)
      external <- data.frame(
        jaccard_index          = jaccard_index(clusters, labels),
        adjusted_rand_index    = adjusted_rand_index(clusters, labels),
        normalized_mutual_info = normalized_mutual_info(clusters, labels)
      )

      metrics_table <- rbind(metrics_table, data.frame(
        dataset     = dataset_name,
        algorithm   = algorithm_name,
        n_clusters  = length(unique(clusters)),
        noise_pct   = round(100 * mean(clusters == 0), 2),
        cbind(internal, external)
      ))

      # vizualizacija (2D ašys – rinkinio; kitu atveju – pirmos dvi PK)
      visual <- plot_cluster_visual(dataset_name, algorithm_name, viz, clusters,
                                    axis_is_pc = (ncol > 2))
      plot_index <- plot_index + 1
      all_plots[[plot_index]] <- visual
      ggsave(file.path(settings_cluster_all$visual_dir,
                       sprintf("%02d_%s_%s.png",
                               plot_index, algorithm_name,
                               sub("\\.csv$", "", dataset_name))),
             visual,
             width  = settings_cluster_all$plot$width,
             height = settings_cluster_all$plot$height,
             dpi    = settings_cluster_all$plot$dpi, units = "in")
    }
  }

  # bendra 3x5 grafikų matrica
  matrix_plot <- arrangeGrob(grobs = all_plots, ncol = 3)
  ggsave(file.path(settings_cluster_all$visual_dir, "clustering_matrix.png"),
         matrix_plot,
         width  = settings_cluster_all$plot$matrix_width,
         height = settings_cluster_all$plot$matrix_height,
         dpi    = settings_cluster_all$plot$dpi, units = "in")

  # matų lentelė
  metrics_path <- file.path(settings_cluster_all$visual_dir,
                            settings_cluster_all$metrics_file)
  write.csv(metrics_table, metrics_path, row.names = FALSE)

  say(verbose, "\n  Matų lentelė:\n")
  print(metrics_table, digits = 3, row.names = FALSE)
  say(verbose, "\n  Grafikai:", settings_cluster_all$visual_dir, "\n")
  say(verbose, "  Matų CSV:", metrics_path, "\n")

  invisible(list(metrics = metrics_table, plots = all_plots))
}