library(ggplot2)

n_runs <- 20L
sample_fraction <- 0.80
k <- 3L
nstart <- 25L
iter_max <- 100L
kmeans_seed <- 29L

umap_config <- umap::umap.defaults
umap_config$n_neighbors <- 15L
umap_config$min_dist <- 0.1
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200L
umap_config$random_state <- 42L

data <- read.csv("a05_unscaled.csv")
features <- data[, names(data) != "class"]
spaces <- c("Original", "PCA 2D", "PCA 6D", "UMAP 2D", "UMAP 6D")

assignments <- array(NA_integer_, c(nrow(features), n_runs + 1L, length(spaces)),
                     dimnames = list(NULL, NULL, spaces))
internal <- array(NA_real_, c(n_runs + 1L, 3L, length(spaces)),
                  dimnames = list(NULL, c("Silhouette", "CH", "DB"), spaces))

cat("Date:", format(Sys.Date()), "| Scaling: z-score within each fit\n")
cat(sprintf("%d samples of %.0f%%, sampling seeds 0-%d\n",
            n_runs, 100 * sample_fraction, n_runs - 1L))
cat(sprintf("Hartigan-Wong: k=%d, nstart=%d, iter.max=%d, seed=%d\n",
            k, nstart, iter_max, kmeans_seed))
cat(sprintf("UMAP: neighbors=%d, min_dist=%.1f, metric=%s, epochs=%d, seed=%d\n",
            umap_config$n_neighbors, umap_config$min_dist, umap_config$metric,
            umap_config$n_epochs, umap_config$random_state))

# Run 0 is the full-data reference. Runs 1-20 use shared subsamples.
for (run in 0:n_runs) {
  if (run == 0L) {
    ids <- seq_len(nrow(features))
    message("Full-data reference")
  } else {
    set.seed(run - 1L)
    ids <- sort(sample.int(nrow(features), round(sample_fraction * nrow(features))))
    message(sprintf("Sample %d/%d", run, n_runs))
  }

  scaled <- scale(features[ids, , drop = FALSE])
  pca <- prcomp(scaled, rank. = 6L)$x
  umap_config$n_components <- 2L
  umap2d <- umap::umap(scaled, config = umap_config)$layout
  umap_config$n_components <- 6L
  umap6d <- umap::umap(scaled, config = umap_config)$layout
  datasets <- list(scaled, pca[, 1:2, drop = FALSE], pca, umap2d, umap6d)

  for (i in seq_along(spaces)) {
    x <- datasets[[i]]
    set.seed(kmeans_seed)
    fit <- kmeans(x, centers = k, nstart = nstart, iter.max = iter_max,
                  algorithm = "Hartigan-Wong")

    assignments[ids, run + 1L, i] <- fit$cluster
    sil <- cluster::silhouette(fit$cluster, dist(x))
    criteria <- clusterCrit::intCriteria(x, fit$cluster,
                                         c("Calinski_Harabasz", "Davies_Bouldin"))
    internal[run + 1L, , i] <- c(summary(sil)$avg.width,
                                 criteria$calinski_harabasz, criteria$davies_bouldin)
  }
}

partition_scores <- function(a, b) {
  c(ARI = aricode::ARI(a, b),
    NMI = aricode::NMI(a, b),
    Jaccard = clusterCrit::extCriteria(a, b, "Jaccard")$jaccard)
}

# Compare the cluster labels of diamonds present in both samples.
pairs <- combn(seq_len(n_runs), 2L)
pair_scores <- data.frame(
  Space = rep(spaces, each = ncol(pairs)),
  First = rep(pairs[1L, ], length(spaces)),
  Second = rep(pairs[2L, ], length(spaces)),
  ARI = NA_real_, NMI = NA_real_, Jaccard = NA_real_
)
for (i in seq_len(nrow(pair_scores))) {
  a <- assignments[, pair_scores$First[i] + 1L, pair_scores$Space[i]]
  b <- assignments[, pair_scores$Second[i] + 1L, pair_scores$Space[i]]
  common <- !is.na(a) & !is.na(b)
  pair_scores[i, c("ARI", "NMI", "Jaccard")] <- partition_scores(a[common], b[common])
}
pair_scores$Space <- factor(pair_scores$Space, levels = spaces)

results <- data.frame(Space = spaces, ARI = "", NMI = "", Jaccard = "",
                      Silhouette = "", CH = "", DB = "")
full_results <- data.frame(Space = spaces, t(internal[1L, , ]),
                           ARI = NA_real_, NMI = NA_real_, Jaccard = NA_real_)
for (i in seq_along(spaces)) {
  pairs_in_space <- pair_scores[pair_scores$Space == spaces[i], c("ARI", "NMI", "Jaccard")]
  scores <- c(
    sprintf("%.3f +/- %.3f", colMeans(pairs_in_space), vapply(pairs_in_space, sd, numeric(1))),
    sprintf("%.3f +/- %.3f", colMeans(internal[-1L, , i]), apply(internal[-1L, , i], 2L, sd))
  )
  results[i, -1L] <- scores
  full_results[i, c("ARI", "NMI", "Jaccard")] <-
    partition_scores(assignments[, 1L, "Original"], assignments[, 1L, i])
}

cat("\nStability and internal metrics (mean +/- SD)\n")
print(results, row.names = FALSE)
cat("\nFull-data fits and agreement with the original space\n")
print(full_results, row.names = FALSE, digits = 4)

convergence <- data.frame()
for (n in c(5L, 10L, 15L, 20L)) {
  means <- aggregate(ARI ~ Space, pair_scores[pair_scores$Second <= n, ], mean)
  convergence <- rbind(convergence, transform(means, Runs = n))
}

print(ggplot(pair_scores, aes(Space, ARI)) +
        geom_boxplot() + labs(x = NULL, title = "ARI across sample pairs") +
        theme_minimal())

print(ggplot(convergence, aes(Runs, ARI, colour = Space)) +
        geom_line() + geom_point() +
        scale_x_continuous(breaks = c(5L, 10L, 15L, 20L)) +
        labs(x = "Number of samples", y = "Mean ARI", title = "Stability convergence") +
        theme_minimal())
