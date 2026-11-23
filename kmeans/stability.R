n_runs <- 20L
sample_fraction <- 0.80

data <- read.csv("a05_unscaled.csv")
features <- data[, names(data) != "class"]
sample_size <- round(sample_fraction * nrow(features))

methods <- data.frame(
  Method = c("PCA", "PCA", "UMAP", "UMAP"),
  Dimensions = c(2L, 6L, 2L, 6L)
)
metrics <- c("ARI", "NMI", "Jaccard", "Silhouette", "CH", "DB")

umap_config <- umap::umap.defaults
umap_config$n_neighbors <- 15L
umap_config$min_dist <- 0.1
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200L
umap_config$random_state <- 42L

reduce_sample <- function(scaled, method, dimensions) {
  if (method == "PCA") {
    return(prcomp(scaled, rank. = dimensions)$x)
  }

  config <- umap_config
  config$n_components <- dimensions
  umap::umap(scaled, config = config)$layout
}

cluster_sample <- function(x) {
  set.seed(29L)
  kmeans(x, centers = 3L, nstart = 25L, iter.max = 100L)$cluster
}

summaries <- vector("list", nrow(methods))

for (i in seq_len(nrow(methods))) {
  method <- methods$Method[i]
  dimensions <- methods$Dimensions[i]
  space <- paste(method, paste0(dimensions, "D"))
  scores <- matrix(NA_real_, nrow = n_runs, ncol = length(metrics),
                   dimnames = list(NULL, metrics))

  for (run in seq_len(n_runs)) {
    message(sprintf("%s, run %d/%d", space, run, n_runs))
    set.seed((run - 1L) * nrow(methods) + i)
    rows_a <- sample.int(nrow(features), sample_size)
    rows_b <- sample.int(nrow(features), sample_size)

    a <- scale(features[rows_a, ])
    b <- scale(features[rows_b, ])
    b <- reduce_sample(b, method, dimensions)
    clusters_a <- cluster_sample(a)
    clusters_b <- cluster_sample(b)

    shared_rows <- intersect(rows_a, rows_b)
    shared_a <- clusters_a[match(shared_rows, rows_a)]
    shared_b <- clusters_b[match(shared_rows, rows_b)]
    scores[run, "ARI"] <- aricode::ARI(shared_a, shared_b)
    scores[run, "NMI"] <- aricode::NMI(shared_a, shared_b)
    scores[run, "Jaccard"] <- clusterCrit::extCriteria(shared_a, shared_b, "Jaccard")$jaccard

    silhouette <- cluster::silhouette(clusters_b, dist(b))
    internal <- clusterCrit::intCriteria(
      b, clusters_b, c("Calinski_Harabasz", "Davies_Bouldin")
    )
    scores[run, "Silhouette"] <- summary(silhouette)$avg.width
    scores[run, "CH"] <- internal$calinski_harabasz
    scores[run, "DB"] <- internal$davies_bouldin
  }

  summaries[[i]] <- data.frame(
    Space = space,
    Metric = metrics,
    Mean = colMeans(scores),
    SD = apply(scores, 2L, sd),
    row.names = NULL
  )
}

results <- do.call(rbind, summaries)
write.csv(results, "stability_results.csv", row.names = FALSE)
print(results, row.names = FALSE)
