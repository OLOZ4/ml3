# =============================================================================
#  UŽDUOTIS 5: stabilumo eksperimentas — imčių pakartojimų analizė
#
#  Kiekvienoje imtyje kartojamas VISAS apdorojimas (ne tik klasterizavimas!):
#    1. 80 % atsitiktinė imtis (stulpelių eilutės)
#    2. z-skaliavimas (RobustScaler analogas — kadangi visas konvejeris
#       jau z-skaliuoja, imties viduje skaliuojama atskirai)
#    3. Dimensijų mažinimas (PCA / UMAP) — kiekvienai imčiai IŠ NAUJO,
#       kitaip visos imtys paveldėtų tą pačią projekciją ir stabilumas
#       būtų pervertintas
#    4. Klasterizavimas (k-means / hierarchinis / DBSCAN)
#
#  Stabilumas: visų imčių porų ARI/NMI/Jaccard TIK BENDRIEMS objektams
#  (sankirtos būdas — tinka visiems algoritmams, nes hierarchinis ir
#  DBSCAN neturi predict). Objektai derinami pagal ID (source eilutę),
#  ne pagal eilės numerį imtyje. Dvi 80 % imtys vidutiniškai turi
#  ~64 % bendrų objektų.
#
#  Pakankamumas: braižomas ARI vidurkis didėjant N (konvergencijos grafikas);
#  lentelėse — vidurkis ± st. nuokrypis; boxplot ARI pasiskirstymui.
#
#  random_state = 0..N-1: skirtingi, bet užfiksuoti — atkartojama.
#
#  Pagrindinė funkcija: run_stability_experiment()
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cluster)
  library(dbscan)
  library(FNN)
})

if (!exists("say")) say <- function(verbose, ...) if (verbose) cat(...)

# ---- šio failo nustatymai (keisk čia) ----------------------------------------
settings_stability <- list(
  verbose = TRUE,

  n_samples = 10,                # imčių skaičius (N; tikrinkite 5->10->15->20)
  sample_fraction = 0.80,        # imties dalis (80 %)

  # dimensijų mažinimas kiekvienoje imtyje: pca / umap / abu
  reduction = list(methods = c("pca", "umap"),
                    dims = c(2, 6),
                    # umap epochos mažinamos stabilumo eksperimentui
                    # (numatytasis umap.defaults n_epochs ~ 500 -> lėta)
                    umap_epochs = 100),

  # klasterizavimo algoritmų parametrai (tie patys kaip UŽDUOTYJE 4)
  cluster_count = 4,
  kmeans = list(nstart = 25),
  hierarchical = list(method = "ward.D2"),
  dbscan = list(minPts = 50, eps = 1.73),

  # pakankamumo kreivei tikrinami imčių taškai
  convergence_points = c(5, 10, 15),

  data_dir = "datasets",
  visual_dir = "visuals/stability",
  plot = list(width = 10, height = 8, dpi = 150)
)

# ==============================================================================
# POSISTEMĖS: apdorojimo žingsnis VIENAI imčiai (atkarpa iš UŽDUOTIES 4)
# ==============================================================================

# z-skaliavimas (RobustScaler atitikmuo: mediana/IQR; čia z-skaliavimas,
# nes visas konvejeris z-skaliuoja)
scale_sample <- function(x) {
  scale(as.matrix(x))
}

# dimensijų mažinimas imčiai (kviečiamas KIEKVIENAI imčiai atskirai!)
reduce_sample <- function(x_scaled, method, dims) {
  if (method == "pca") {
    prcomp(x_scaled, rank. = dims)$x[, 1:dims, drop = FALSE]
  } else if (method == "umap") {
    cfg <- umap::umap.defaults
    cfg$n_components <- dims
    cfg$seed <- 42
    cfg$n_epochs <- settings_stability$reduction$umap_epochs
    as.matrix(umap::umap(x_scaled, cfg)$layout)
  } else {
    x_scaled
  }
}

# klasterizavimas imčiai (ta patys parametrai kaip UŽDUOTYJE 4)
cluster_sample <- function(x, algorithm) {
  s <- settings_stability
  if (algorithm == "kmeans") {
    set.seed(75)
    kmeans(x, centers = s$cluster_count,
           nstart = s$kmeans$nstart, iter.max = 100)$cluster
  } else if (algorithm == "hierarchical") {
    cutree(hclust(dist(x), method = s$hierarchical$method),
           k = 2)
  } else if (algorithm == "dbscan") {
    as.integer(dbscan(x, eps = s$dbscan$eps,
                       minPts = s$dbscan$minPts)$cluster)
  }
}

# ==============================================================================
# POSISTEMĖS: išorinių matų skaičiavimas TIK BENDRIEMS objektams
# ==============================================================================

# bendrų objektų žymių poros: deriname pagal pirminį indeksą (obj_id)
labels_common <- function(labels_a, ids_a, labels_b, ids_b) {
  common <- intersect(ids_a, ids_b)
  list(a = labels_a[match(common, ids_a)],
       b = labels_b[match(common, ids_b)],
       common_ids = common)
}

# ARI/NMI/Jaccard porai (naudojamos UŽDUOTIES 4 formulės, jei jau įkeltos;
# čia savarankiškos kopijos, kad užduotis būtų savarankiška)
if (!exists("pair_counts")) {
  pair_counts <- function(clusters, labels) {
    contingency <- table(clusters, labels)
    n <- sum(contingency)
    a <- sum(choose(contingency, 2))
    b <- sum(choose(rowSums(contingency), 2)) - a
    c <- sum(choose(colSums(contingency), 2)) - a
    total <- choose(n, 2)
    list(a = a, b = b, c = c, total = total)
  }
}

pairwise_metrics <- function(labels_a, labels_b) {
  counts <- pair_counts(labels_a, labels_b)
  jaccard <- counts$a / (counts$a + counts$b + counts$c)
  expected <- (counts$a + counts$b) * (counts$a + counts$c) / counts$total
  max_index <- (counts$a + counts$b + counts$a + counts$c) / 2
  ari <- if (max_index == expected) NA_real_ else
    (counts$a - expected) / (max_index - expected)
  # NMI
  contingency <- table(labels_a, labels_b)
  n <- sum(contingency)
  p_ij <- contingency / n
  p_i  <- rowSums(contingency) / n
  p_j  <- colSums(contingency) / n
  mi <- sum(p_ij * log(ifelse(p_ij > 0, p_ij, 1) / (p_i %*% t(p_j))))
  h_i <- -sum(p_i * log(ifelse(p_i > 0, p_i, 1)))
  h_j <- -sum(p_j * log(ifelse(p_j > 0, p_j, 1)))
  nmi <- if (h_i * h_j == 0) NA_real_ else 2 * mi / (h_i + h_j)
  c(ari = ari, nmi = nmi, jaccard = jaccard)
}

# ==============================================================================
# PAGRINDINĖ UŽDUOTIES 5 FUNKCIJA
# ==============================================================================

run_stability_experiment <- function() {

  s <- settings_stability
  verbose <- s$verbose
  say(verbose, "=== UŽDUOTIS 5: stabilumo eksperimentas ===\n")
  say(verbose, sprintf("  imčių N = %d, po %.0f%%, algoritmai: k-means/hierarchinis/DBSCAN\n",
                       s$n_samples, 100 * s$sample_fraction))

  dir.create(s$visual_dir, showWarnings = FALSE, recursive = TRUE)

  # bazinis rinkinys (originali erdvė)
  base <- read.csv(file.path(s$data_dir, "a05.csv"), check.names = FALSE)
  base_features <- as.matrix(base[, setdiff(names(base), "class")])
  n <- nrow(base_features)
  algorithms <- c("kmeans", "hierarchical", "dbscan")

  # erdvės: originali + visos reduction kombinacijos
  spaces <- list(c(method = "none", dims = 0))
  for (m in s$reduction$methods)
    for (d in s$reduction$dims)
      spaces <- c(spaces, list(c(method = m, dims = d)))

  space_names <- vapply(spaces, function(sp)
    if (sp["method"] == "none") "original" else
      sprintf("%s_%d", sp["method"], as.integer(sp["dims"])), "")

  # imtys: kiekvienai random_state = 0..N-1
  sample_ids <- lapply(0:(s$n_samples - 1), function(rs) {
    set.seed(rs)
    sort(sample(n, round(n * s$sample_fraction)))
  })

  # rezultatų struktūros
  pair_metrics <- data.frame()   # visos imčių poros x erdvės x algoritmai
  internal_by_sample <- data.frame()  # vidiniai matai kiekvienai imčiai

  # skaičiavimas: erdvė x algoritmas -> N imtys -> visos poros
  for (sp_idx in seq_along(spaces)) {

    sp <- spaces[[sp_idx]]
    sp_name <- space_names[sp_idx]
    say(verbose, sprintf("\n  Erdvė: %s\n", sp_name))

    for (algorithm in algorithms) {

      say(verbose, sprintf("    algoritmas: %s ...\n", algorithm))

      # klasterizuoja KIEKVIENĄ imtį atskirai
      sample_results <- lapply(seq_along(sample_ids), function(i) {
        ids <- sample_ids[[i]]
        x <- base_features[ids, , drop = FALSE]
        x <- scale_sample(x)                 # skalė imties viduje
        if (sp["method"] != "none")
          x <- reduce_sample(x, sp["method"], as.integer(sp["dims"]))
        clusters <- cluster_sample(x, algorithm)
        list(ids = ids, clusters = clusters, x = x)
      })

      # vidiniai matai kiekvienai imčiai
      for (i in seq_along(sample_results)) {
        res <- sample_results[[i]]
        sil <- if (length(unique(res$clusters)) < 2) NA_real_ else
          mean(cluster::silhouette(res$clusters, dist(res$x))[, "sil_width"])
        internal_by_sample <- rbind(internal_by_sample, data.frame(
          space = sp_name, algorithm = algorithm, sample = i - 1,
          silhouette = sil
        ))
      }

      # visos imčių poros: matai TIK bendriems objektams
      for (i in 1:(s$n_samples - 1)) for (j in (i + 1):s$n_samples) {
        ra <- sample_results[[i]]; rb <- sample_results[[j]]
        common <- labels_common(ra$clusters, ra$ids,
                                rb$clusters, rb$ids)
        if (length(common$common_ids) < 2) next
        m <- pairwise_metrics(common$a, common$b)
        pair_metrics <- rbind(pair_metrics, data.frame(
          space = sp_name, algorithm = algorithm,
          pair_i = i - 1, pair_j = j - 1,
          n_common = length(common$common_ids),
          ari = m["ari"], nmi = m["nmi"], jaccard = m["jaccard"]
        ))
      }
    }
  }

  # ==============================================================================
  # IŠVESTIS: lentelės (vidurkis ± SN), boxplot, konvergencijos kreivė
  # ==============================================================================

  # 1. lentelė: vidurkis ± st. nuokrypis visoms kombinacijoms
  summary_table <- do.call(rbind, lapply(split(pair_metrics,
    paste(pair_metrics$space, pair_metrics$algorithm)), function(df) {
    data.frame(
      space = df$space[1], algorithm = df$algorithm[1],
      n_pairs = nrow(df),
      ari_mean = mean(df$ari),
      ari_sd = sd(df$ari),
      nmi_mean = mean(df$nmi),
      nmi_sd = sd(df$nmi),
      jaccard_mean = mean(df$jaccard),
      jaccard_sd = sd(df$jaccard)
    )
  }))

  write.csv(summary_table,
            file.path(s$visual_dir, "stability_summary.csv"),
            row.names = FALSE)

  # 2. boxplot: ARI pasiskirstymas visoms kombinacijoms
  pair_metrics$space <- factor(pair_metrics$space, levels = space_names)
  boxplot_ari <- ggplot(pair_metrics,
                        aes(x = space, y = ari, fill = algorithm)) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4) +
    scale_fill_brewer(palette = "Set2", name = "Algoritmas") +
    labs(title = "Stabilumas: ARI pasiskirstymas (imčių poros)",
         subtitle = sprintf("%d imtys po %.0f%% -> %d poros; tik bendrieji objektai",
                            s$n_samples, 100 * s$sample_fraction,
                            s$n_samples * (s$n_samples - 1) / 2),
         x = "Erdvė", y = "ARI") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(s$visual_dir, "stability_ari_boxplot.png"), boxplot_ari,
         width = s$plot$width, height = s$plot$height,
         dpi = s$plot$dpi, units = "in")

  # 3. konvergencijos kreivė: ARI vidurkis didėjant N (visi algoritmai kartu)
  conv_by_algo <- data.frame()
  for (N in s$convergence_points) {
    ids_keep <- 0:(N - 1)
    df_sub <- pair_metrics[pair_metrics$pair_i %in% ids_keep &
                           pair_metrics$pair_j %in% ids_keep, ]
    for (algo in algorithms) {
      for (sp_name in unique(df_sub$space)) {
        df_a <- df_sub[df_sub$algorithm == algo & df_sub$space == sp_name, ]
        if (!nrow(df_a)) next
        conv_by_algo <- rbind(conv_by_algo, data.frame(
          N = N,
          space = sp_name,
          algorithm = algo,
          ari_mean = mean(df_a$ari),
          ari_sd = sd(df_a$ari)
        ))
      }
    }
  }

  conv_plot <- ggplot(conv_by_algo,
                      aes(x = N, y = ari_mean, colour = algorithm)) +
    geom_point(size = 2.5) +
    geom_line(linewidth = 0.7) +
    geom_errorbar(aes(ymin = ari_mean - ari_sd,
                      ymax = ari_mean + ari_sd), width = 0.4, alpha = 0.5) +
    facet_wrap(~ space) +
    scale_colour_brewer(palette = "Set2", name = "Algoritmas") +
    labs(title = "Pakankamumas: ARI vidurkis didėjant imčių skaičių N",
         subtitle = "Kai vidurkis nusistovi, imčių pakanka",
         x = "Imčių skaičius N", y = "ARI vidurkis") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(s$visual_dir, "stability_convergence.png"), conv_plot,
         width = s$plot$width, height = s$plot$height,
         dpi = s$plot$dpi, units = "in")

  # 4. vidiniai matai lentelėje: vidurkis ± SN kiekvienai erdvei
  internal_summary <- do.call(rbind, lapply(split(internal_by_sample,
    paste(internal_by_sample$space, internal_by_sample$algorithm)),
    function(df) data.frame(
      space = df$space[1], algorithm = df$algorithm[1],
      silhouette_mean = mean(df$silhouette),
      silhouette_sd = sd(df$silhouette)
    )))

  write.csv(internal_summary,
            file.path(s$visual_dir, "stability_internal.csv"),
            row.names = FALSE)

  # spausdinimas
  say(verbose, "\n  Stabilumo lentelė (ARI/NMI/Jaccard vidurkis ± st. nuokrypis):\n")
  print(summary_table, row.names = FALSE, digits = 3)
  say(verbose, "\n  Vidinių matų lentelė (Silhouette vidurkis ± SN per imtis):\n")
  print(internal_summary, row.names = FALSE, digits = 3)
  say(verbose, "\n  Grafikai ir lentelės:", s$visual_dir, "\n")

  invisible(list(pairs = pair_metrics, internal = internal_by_sample,
                 summary = summary_table))
}