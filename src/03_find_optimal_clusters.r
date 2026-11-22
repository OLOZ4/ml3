# =============================================================================
#  UŽDUOTIS 3: rasti ir pagrįsti optimalų klasterių skaičių
#
#  Keturi metodai (kiekvienas vizualizuojamas ir rašomas į visuals/<aplankas>/):
#    1. Alkūnės metodas (elbow)      – WSS kreivė; alkūnė = lūžio taškas
#    2. Silueto metodas (silhouette) – vidutinis silueto plotis kiekvienam k;
#                                      optimalus k = didžiausias plotis
#    3. Empirinis metodas            – heuristinė formulė k = sqrt(n / 2)
#    4. Dendrograma                  – hierarchinis klasterizavimas; k = kliauties
#                                      vietoje, kur didžiausias susiliejimo aukščio
#                                      šuolis
#
#  (Kelio (knee) metodas čia nenaudojamas – jis skirtas DBSCAN eps parinkimui,
#   žr. src/04_cluster_datasets.r)
#
#  Pagrindinė funkcija: find_optimal_clusters(data)
#  Posistemės: compute_elbow(), compute_silhouette(),
#              compute_empirical(), compute_dendrogram(),
#              plot_elbow(), plot_silhouette(), plot_empirical(),
#              plot_dendrogram(), plot_comparison(), save_plot()
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cluster))

# pagalbinė: spausdinti tik kai verbose = TRUE
if (!exists("say")) say <- function(verbose, ...) if (verbose) cat(...)

# ---- šio failo nustatymai (keisk čia) ----------------------------------------
settings_clusters <- list(
  verbose = TRUE,                 # TRUE = spausdinti šios užduoties žinutes

  k_min = 2,                      # mažiausias tikrinamas klasterių skaičius
  k_max = 10,                     # didžiausias tikrinamas klasterių skaičius
  nstart = 25,                    # k-means paleidimų su atsitiktiniais centrais
  seed = 75,                      # atkartojamumas

  # alkūnės aptikimo taisyklė: mažiausias k, kurio WSS pagerėjimas
  # ( santykinis ) jau mažesnis už šią ribą
  elbow_drop_threshold = 0.10,

  # dendrogramai: hierarchinio klasterizavimo susiejimo metodas
  hclust_method = "ward.D2",

  visual_dir = "visuals/optimal_clusters",   # kur rašyti grafikus
  plot = list(width = 10, height = 8, dpi = 150)   # visi grafikai vienodi
)

# ------------------------------------------------------------------------------
# posistemė: k-means pagal nustatymus
run_kmeans <- function(data, k) {
  set.seed(settings_clusters$seed)
  kmeans(data, centers = k, nstart = settings_clusters$nstart,
         iter.max = 100)
}

# posistemė: alkūnės metodas – WSS kiekvienam k + alkūnės taškas
compute_elbow <- function(features) {
  k_values <- settings_clusters$k_min : settings_clusters$k_max
  wss <- sapply(k_values, function(k) run_kmeans(features, k)$tot.withinss)
  result <- data.frame(k = k_values, wss = wss)

  # alkūnė = pirmas k, kuriuo santykinis WSS sumažėjimas < ribos
  relative_drop <- c(NA, -diff(result$wss) / head(result$wss, -1))
  first_below <- which(relative_drop < settings_clusters$elbow_drop_threshold)[1]
  optimal_k <- if (is.na(first_below)) result$k[which.min(result$wss)]
               else result$k[first_below - 1]

  list(curve = result, relative_drop = relative_drop, optimal_k = optimal_k)
}

# posistemė: kelio (knee) metodas – WSS kreivės taškas, toliausiai nuo
# tiesės, jungiančios kreivės pradžią ir galą (naudoja tą pačią WSS kreivę)
compute_knee <- function(elbow) {
  curve <- elbow$curve
  # tiesė nuo pirmo iki paskutinio kreivės taško
  slope <- (curve$wss[nrow(curve)] - curve$wss[1]) /
           (curve$k[nrow(curve)] - curve$k[1])
  line_y <- curve$wss[1] + slope * (curve$k - curve$k[1])
  distances <- abs(curve$wss - line_y)
  knee_k <- curve$k[which.max(distances)]
  list(curve = curve, line_y = line_y, distances = distances,
       optimal_k = knee_k)
}

plot_knee <- function(knee) {
  curve <- knee$curve
  optimal <- knee$optimal_k
  line_data <- data.frame(k = curve$k, line = knee$line_y)
  ggplot(curve, aes(x = k, y = wss)) +
    geom_line(colour = "#2E86AB", linewidth = 0.7) +
    geom_point(size = 2) +
    geom_line(data = line_data, aes(x = k, y = line),
              linetype = "dashed", colour = "grey60", linewidth = 0.5) +
    geom_vline(xintercept = optimal, linetype = "dashed", colour = "#A23B72") +
    annotate("text", x = optimal, y = max(curve$wss),
             label = paste0("kelias: k = ", optimal),
             vjust = -0.5, hjust = -0.1, colour = "#A23B72", fontface = "bold") +
    labs(title = "Kelio metodas (knee)",
         subtitle = "Taškas, toliausiai nuo tiesės tarp kreivės galų",
         x = "Klasterių skaičius k", y = "WSS") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

# posistemė: dendrograma – hierarchinis klasterizavimas; optimalus k =
# klasterių skaičius, kai dendrogramą kirto didžiausiame susiliejimo aukščio šuolyje
compute_dendrogram <- function(features) {
  hierarchy <- hclust(dist(features), method = settings_clusters$hclust_method)
  merge_heights <- hierarchy$height
  # didžiausias šuolis tarp gretimų susiliejimų; kirpdami jame gauname
  # n - j klasterių (j = šuolio indeksas)
  gap_index <- which.max(diff(merge_heights))
  optimal_k <- nrow(features) - gap_index
  list(hierarchy = hierarchy, optimal_k = optimal_k)
}

# posistemė: silueto metodas – vidutinis silueto plotis kiekvienam k
compute_silhouette <- function(features) {
  k_values <- settings_clusters$k_min : settings_clusters$k_max
  distances <- dist(features)
  widths <- sapply(k_values, function(k) {
    clustering <- run_kmeans(features, k)
    mean(cluster::silhouette(clustering$cluster, distances)[, "sil_width"])
  })
  result <- data.frame(k = k_values, avg_width = widths)
  list(curve = result, optimal_k = result$k[which.max(result$avg_width)])
}

# posistemė: empirinis metodas – heuristinė formulė k = sqrt(n / 2)
compute_empirical <- function(n) {
  round(sqrt(n / 2))
}

# ------------------------------------------------------------------------------
# posistemės: grafikai

save_plot <- function(plot, file_name) {
  file_path <- file.path(settings_clusters$visual_dir, file_name)
  ggsave(file_path, plot,
         width  = settings_clusters$plot$width,
         height = settings_clusters$plot$height,
         dpi    = settings_clusters$plot$dpi, units = "in")
  say(settings_clusters$verbose, "  grafikas:", file_path, "\n")
}

plot_elbow <- function(elbow) {
  optimal <- elbow$optimal_k
  ggplot(elbow$curve, aes(x = k, y = wss)) +
    geom_line(colour = "#2E86AB", linewidth = 0.7) +
    geom_point(size = 2) +
    geom_vline(xintercept = optimal, linetype = "dashed", colour = "#A23B72") +
    annotate("text", x = optimal, y = max(elbow$curve$wss),
             label = paste0("alkūnė: k = ", optimal),
             vjust = -0.5, hjust = -0.1, colour = "#A23B72", fontface = "bold") +
    labs(title = "Alkūnės metodas (elbow)",
         subtitle = "Klasterių vidinė kvadratų suma (WSS) pagal klasterių skaičių",
         x = "Klasterių skaičius k", y = "WSS") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

plot_silhouette <- function(sil) {
  optimal <- sil$optimal_k
  ggplot(sil$curve, aes(x = k, y = avg_width)) +
    geom_col(fill = ifelse(sil$curve$k == optimal, "#A23B72", "#2E86AB")) +
    labs(title = "Silueto metodas",
         subtitle = paste0("Vidutinis silueto plotis; optimalus k = ", optimal),
         x = "Klasterių skaičius k", y = "Vidutinis silueto plotis") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

# dendrograma piešiama baziniais grafikais (ne ggplot) – naudojame png()
plot_dendrogram <- function(dendro, file_name) {
  file_path <- file.path(settings_clusters$visual_dir, file_name)
  p <- settings_clusters$plot
  png(file_path, width = p$width, height = p$height,
      units = "in", res = p$dpi)
  plot(dendro$hierarchy, labels = FALSE,
       main = "Dendrograma (hierarchinis klasterizavimas)",
       xlab = "Objektai", ylab = "Susiliejimo aukštis")
  rect.hclust(dendro$hierarchy, k = dendro$optimal_k,
              border = "#A23B72")
  dev.off()
  say(settings_clusters$verbose, "  grafikas:", file_path, "\n")
  invisible(NULL)
}

plot_empirical <- function(n) {
  # formulės kreivė: k = sqrt(n/2) priklausomai nuo duomenų dydžio
  sizes <- seq(100, 10000, by = 100)
  curve_data <- data.frame(n = sizes, k = sqrt(sizes / 2))
  our_k <- compute_empirical(n)
  ggplot(curve_data, aes(x = n, y = k)) +
    geom_line(colour = "#2E86AB", linewidth = 0.7) +
    geom_point(data = data.frame(n = n, k = our_k),
               aes(x = n, y = k), size = 3, colour = "#A23B72") +
    annotate("text", x = n, y = our_k,
             label = paste0("mūsų duomenys: n = ", n, " -> k = ", our_k),
             hjust = 1.05, vjust = -1, colour = "#A23B72", fontface = "bold") +
    labs(title = "Empirinis metodas",
         subtitle = "Heuristinė formulė k = sqrt(n / 2)",
         x = "Objektų skaičius n", y = "Empirinis klasterių skaičius k") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

plot_comparison <- function(estimates) {
  ggplot(estimates, aes(x = method, y = k, fill = method)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = k), vjust = -0.5, fontface = "bold") +
    labs(title = "Metodų palyginimas",
         subtitle = "Kiekvieno metodo siūlomas klasterių skaičius",
         x = NULL, y = "Siūlomas klasterių skaičius k") +
    scale_fill_manual(values = c("Alkūnė" = "#2E86AB",
                                 "Kelias" = "#5B8FF9",
                                 "Siluetas" = "#A23B72",
                                 "Empirinis" = "#F18F01",
                                 "Dendrograma" = "#7C9E4F")) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

# ------------------------------------------------------------------------------
# pagrindinė užduoties 3 funkcija
find_optimal_clusters <- function(data) {

  verbose <- settings_clusters$verbose

  feature_names <- setdiff(names(data), "class")
  features <- data[, feature_names, drop = FALSE]

  dir.create(settings_clusters$visual_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. alkūnės metodas
  elbow <- compute_elbow(features)
  say(verbose, sprintf("  alkūnės metodas: optimalus k = %d\n", elbow$optimal_k))
  save_plot(plot_elbow(elbow), "01_elbow.png")

  # 2. kelio (knee) metodas – naudoja tą pačią WSS kreivę
  knee <- compute_knee(elbow)
  say(verbose, sprintf("  kelio metodas: optimalus k = %d\n", knee$optimal_k))
  save_plot(plot_knee(knee), "02_knee.png")

  # 3. silueto metodas
  sil <- compute_silhouette(features)
  say(verbose, sprintf("  silueto metodas: optimalus k = %d (plotis %.4f)\n",
                       sil$optimal_k,
                       max(sil$curve$avg_width)))
  save_plot(plot_silhouette(sil), "03_silhouette.png")

  # 4. empirinis metodas
  empirical_k <- compute_empirical(nrow(features))
  say(verbose, sprintf("  empirinis metodas (sqrt(n/2)): k = %d\n", empirical_k))
  save_plot(plot_empirical(nrow(features)), "04_empirical.png")

  # 5. dendrograma (hierarchinis klasterizavimas)
  say(verbose, "  dendrograma: skaičiuojamas hierarchinis klasterizavimas...\n")
  dendro <- compute_dendrogram(features)
  say(verbose, sprintf("  dendrograma: optimalus k = %d (didžiausias aukščio šuolis)\n",
                       dendro$optimal_k))
  plot_dendrogram(dendro, "05_dendrogram.png")

  # metodų palyginimas viename grafike
  estimates <- data.frame(
    method = c("Alkūnė", "Kelias", "Siluetas", "Empirinis", "Dendrograma"),
    k      = c(elbow$optimal_k, knee$optimal_k, sil$optimal_k,
               empirical_k, dendro$optimal_k)
  )
  save_plot(plot_comparison(estimates), "06_comparison.png")

  say(verbose, "  metodų įverčiai:", paste(estimates$method, estimates$k,
                                           collapse = " | "), "\n")

  invisible(list(elbow = elbow, knee = knee, silhouette = sil,
                 empirical_k = empirical_k, dendrogram = dendro,
                 estimates = estimates))
}