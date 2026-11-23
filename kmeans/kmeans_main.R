library(factoextra)

elbow <- function(data, nclust) {
  set.seed(6)
  wcss = vector()
  for (i in 1:nclust) wcss[i] = sum(kmeans(data, i, iter.max = 100)$withinss)
  plot(1:nclust,
       wcss,
       type = 'b',
       main = paste(' Elbow metodas'),
       xlab = 'Klasterių skaičius',
       ylab = 'WCSS')
}

cluster <- function(seed, data, centers, algorithm = "Hartigan-Wong") {
  set.seed(seed)
  kmeans(data, centers = centers, iter.max = 100, nstart = 25,
         algorithm = algorithm)
}

compare_algorithms <- function(data) {
  distances <- dist(data)
  for (algorithm in c("Hartigan-Wong", "Lloyd", "MacQueen")) {
    fit <- cluster(29, data, 3, algorithm)
    sil <- cluster::silhouette(fit$cluster, distances)
    cat(sprintf("%s: WSS %.2f, silhouette %.4f, iterations %d\n",
                algorithm, fit$tot.withinss, summary(sil)$avg.width, fit$iter))
  }
}

# Originali
dataset = read.csv('a05.csv')
dataset <- dataset[1:8]

# PCA 6D
pc6d <- as.matrix(read.csv('a05_pca_6.csv')[1:6])

# PCA 2D
pc2d <- as.matrix(read.csv('a05_pca_2.csv')[1:2])

# UMAP 6D
umap6d <- as.matrix(read.csv('a05_umap_6.csv')[1:6])

# UMAP 2D
umap2d <- as.matrix(read.csv('a05_umap_2.csv')[1:2])

# ELBOW
elbow(dataset, 8)

# SILHOUETTE
fviz_nbclust(dataset, kmeans, method='silhouette')

# KLASTERIZAVIMAS
dataset.res <- cluster(29, dataset, 3)
pc6d.res <- cluster(29, pc6d, 3)
pc2d.res <- cluster(29, pc2d, 3)
umap6d.res <- cluster(29, umap6d, 3)
umap2d.res <- cluster(29, umap2d, 3)

# KLASTERIU POZYMIU VIDURKIAI
# ISVEDAMI IS NESTANDARTIZUOTOS AIBES
dataset_unscaled <- read.csv("a05_unscaled.csv")[names(dataset)]
cluster_means <- aggregate(dataset_unscaled, by = list(Cluster = dataset.res$cluster), FUN = mean)
print(cluster_means, row.names = FALSE)

# SKIRTINGU ALGORITMU PALYGINIMAS
compare_algorithms(dataset)

# KLASTERIU VIZUALIZAVIMAS
library(cluster)
cluster_plots <- list(
  fviz_cluster(dataset.res, dataset, stand = FALSE, labelsize = 0) + ggplot2::labs(title = "Originali aibė"),
  fviz_cluster(pc6d.res, pc6d, choose.vars = c("pca1", "pca2"),
               stand = FALSE, labelsize = 0) + ggplot2::labs(title = "PCA 6D"),
  fviz_cluster(pc2d.res, pc2d, stand = FALSE, labelsize = 0) + ggplot2::labs(title = "PCA 2D"),
  fviz_cluster(umap6d.res, umap6d, stand = FALSE, labelsize = 0) + ggplot2::labs(title = "UMAP 6D"),
  fviz_cluster(umap2d.res, umap2d, stand = FALSE, labelsize = 0) + ggplot2::labs(title = "UMAP 2D")
)

for (cluster_plot in cluster_plots) print(cluster_plot)

combined_cluster_plot <- ggpubr::ggarrange(plotlist = cluster_plots, ncol = 3, nrow = 2)
print(combined_cluster_plot)

# SILUETŲ VIZUALIZAVIMAS
sil <- silhouette(dataset.res$cluster, dist(dataset))
fviz_silhouette(sil) + ggplot2::labs(title = "Originali aibė")

sil_pc6d <- silhouette(pc6d.res$cluster, dist(pc6d))
fviz_silhouette(sil_pc6d) + ggplot2::labs(title = "PCA 6D")

sil_pc2d <- silhouette(pc2d.res$cluster, dist(pc2d))
fviz_silhouette(sil_pc2d) + ggplot2::labs(title = "PCA 2D")

sil_umap6d <- silhouette(umap6d.res$cluster, dist(umap6d))
fviz_silhouette(sil_umap6d) + ggplot2::labs(title = "UMAP 6D")

sil_umap2d <- silhouette(umap2d.res$cluster, dist(umap2d))
fviz_silhouette(sil_umap2d) + ggplot2::labs(title = "UMAP 2D")

# VIDINIAI MATAI
internal_criteria <- function(data, clusters, sil) {
  criteria <- clusterCrit::intCriteria(
    as.matrix(data), clusters, c("Calinski_Harabasz", "Davies_Bouldin")
  )
  criteria$silhouette <- summary(sil)$avg.width
  criteria
}

internal_criteria(dataset, dataset.res$cluster, sil)
internal_criteria(pc6d, pc6d.res$cluster, sil_pc6d)
internal_criteria(pc2d, pc2d.res$cluster, sil_pc2d)
internal_criteria(umap6d, umap6d.res$cluster, sil_umap6d)
internal_criteria(umap2d, umap2d.res$cluster, sil_umap2d)
