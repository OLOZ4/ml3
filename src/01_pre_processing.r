# =============================================================================
#  pre_processing.r — visos apdorojimo funkcijos viename faile
#
#    read_raw_data()               – UŽDUOTIS 1: nuskaitymas, € nuėmimas
#    fill_invalid_values()         – UŽDUOTIS 2: netinkamų reikšmių užpildymas
#    remove_outliers()             – UŽDUOTIS 3: išskirčių šalinimas (IQR)
#    remove_correlated_features() – UŽDUOTIS 4: koreliuotų požymių šalinimas
#    normalize_z_score()           – UŽDUOTIS 5: z-skaliavimas
# =============================================================================

# UŽDUOTIS 1: nuskaito CSV, nuima "€" iš price, tuščias reikšmes paverčia NA.

read_raw_data <- function(input_path) {

  raw <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)

  for (name in setdiff(names(raw), "class")) {
    text <- trimws(as.character(raw[[name]]))
    if (name == "price") text <- trimws(gsub("€", "", text, fixed = TRUE))
    values <- suppressWarnings(as.numeric(text))
    raw[[name]] <- ifelse(is.finite(values), values, NA_real_)
  }

  raw$class <- factor(raw$class)
  cat("Nuskaityta:", nrow(raw), "eilučių,",
      length(setdiff(names(raw), "class")), "požymių\n")
  raw
}

# UŽDUOTIS 2 ------------------------------------------------------------------
# Užpildo trūkstamas / nulines / neigiamas bazes (carat, price, x, y, z):
#   1. pagal kitus rodiklius suskaičiuojami kandidatai (formulės žemiau);
#   2. jei bent 2 kandidatai sutampa (< 1% paklaida) – imama jų mediana;
#   3. kitu atveju: carat/price – tiesinė regresija, x/y/z – klasės mediana;
#   4. išvestiniai požymiai perskaičiuojami iš pataisytų bazių.
# Ratai kartojami iki fiksuoto taško (nauji užpildymai atveria naujų kandidatų).

# kandidatų formulės kiekvienai bazei:
# NA operandas -> NA kandidatas; dalyba iš 0 -> Inf (ne galimas) -> atmetama
candidate_formulas <- list(
  carat = c("price / price_per_carat",
            "carat_per_volume * volume_xyz",
            "price_per_volume * volume_xyz / price_per_carat"),
  price = c("price_per_carat * carat",
            "price_per_volume * volume_xyz",
            "price_per_carat * carat_per_volume * volume_xyz"),
  x     = c("area_xy / y", "area_xz / z", "volume_xyz / (y * z)",
            "length_width_ratio * y", "3 * mean_dimension - y - z"),
  y     = c("area_xy / x", "area_yz / z", "volume_xyz / (x * z)",
            "x / length_width_ratio", "3 * mean_dimension - x - z"),
  z     = c("area_xz / x", "area_yz / y", "volume_xyz / (x * y)",
            "depth_ratio * (x + y) / 2", "3 * mean_dimension - x - y")
)

# ar reikšmės yra teigiami galimi skaičiai?
is_valid <- function(values) is.finite(values) & values > 0

# iš vienos eilutės kandidatų: jei bent 2 sutampa – jų mediana, kitu atveju NA
fill_value_from_candidates <- function(candidates, tolerance = 0.01) {
  valid <- sort(candidates[is_valid(candidates)])
  if (length(valid) < 2) return(NA_real_)
  # paliekame tik tuos, kurie sutampa dar bent su vienu kandidatu
  agreed <- valid[vapply(valid, function(v)
    sum(abs(valid - v) <= tolerance * v), 0L) >= 2]
  if (length(agreed)) median(agreed) else NA_real_
}

# tos pačios klasės galiojančių reikšmių mediana
class_median <- function(data, column, row) {
  values <- data[[column]][data$class == data$class[row]]
  values <- values[is_valid(values)]
  if (length(values)) median(values) else NA_real_
}

# teigiama regresijos prognozė (kitu atveju NA)
predict_positive <- function(model, data, row) {
  prediction <- suppressWarnings(unname(predict(model, newdata = data[row, ])))
  if (is.finite(prediction) && prediction > 0) prediction else NA_real_
}

# išvestinių požymių perskaičiavimas (po užpildymo visos bazės > 0)
recompute_derived_features <- function(data) {
  within(data, {
    volume_xyz         <- x * y * z
    area_xy            <- x * y
    area_xz            <- x * z
    area_yz            <- y * z
    price_per_carat    <- price / carat
    length_width_ratio <- x / y
    depth_ratio        <- z / ((x + y) / 2)
    table_depth_ratio  <- table / depth
    carat_per_volume   <- carat / volume_xyz
    price_per_volume   <- price / volume_xyz
    mean_dimension     <- (x + y + z) / 3
    dimension_cv       <- sqrt(((x - mean_dimension)^2 +
                                 (y - mean_dimension)^2 +
                                 (z - mean_dimension)^2) / 3) / mean_dimension
  })
}

fill_invalid_values <- function(data) {

  base_columns <- c("carat", "price", "x", "y", "z")

  for (round in 1:10) {                       # ratai iki fiksuoto taško

    n_filled <- 0

    # atsarginės regresijos – tik ant galiojančių eilučių
    carat_model <- lm(carat ~ x + y + z,
                      data = data[is_valid(data$carat) & is_valid(data$x) &
                                    is_valid(data$y) & is_valid(data$z), ])
    price_model <- lm(price ~ carat + class,
                      data = data[is_valid(data$price) &
                                    is_valid(data$carat), ])

    fallbacks <- list(
      carat = function(d, r) predict_positive(carat_model, d, r),
      price = function(d, r) predict_positive(price_model, d, r),
      x     = function(d, r) class_median(d, "x", r),
      y     = function(d, r) class_median(d, "y", r),
      z     = function(d, r) class_median(d, "z", r)
    )

    cat("Iteracija", round, ":\n")

    for (column in base_columns) {

      rows_to_fill <- which(!is_valid(data[[column]]))
      if (!length(rows_to_fill)) next

      # visos kandidatų formulės įvertinamos vektoriškai visoms eilutėms
      candidate_table <- sapply(candidate_formulas[[column]],
                                 function(f) eval(parse(text = f), data))
      filled <- 0

      for (r in rows_to_fill) {
        value <- fill_value_from_candidates(as.numeric(candidate_table[r, ]))
        if (is.na(value)) value <- fallbacks[[column]](data, r)
        if (is_valid(value)) {
          data[r, column] <- value
          filled <- filled + 1
        }
      }

      cat(sprintf("  %-6s: užpildyta %d / %d\n",
                  column, filled, length(rows_to_fill)))
      n_filled <- n_filled + filled
    }

    data <- recompute_derived_features(data)
    if (n_filled == 0) break
  }

  # patikra: visos bazės turi būti teigiami skaičiai
  stopifnot(all(vapply(data[base_columns], function(v) all(is_valid(v)), TRUE)))
  data
}

# UŽDUOTIS 3: šalina eilutes, kurių stulpelis viršija Q3 + multiplier * IQR.

remove_outliers <- function(data, column = "price_per_carat", iqr_multiplier = 1.5) {

  q1 <- quantile(data[[column]], 0.25)
  q3 <- quantile(data[[column]], 0.75)
  upper_fence <- q3 + iqr_multiplier * (q3 - q1)

  keep <- data[[column]] <= upper_fence

  cat(sprintf("%s: riba Q3 + %.1fxIQR = %.2f, šaliname %d eilutes (%.1f%%)\n",
              column, iqr_multiplier, upper_fence,
              sum(!keep), 100 * mean(!keep)))
  data[keep, ]
}

# UŽDUOTIS 4: šalina koreliuotus požymius.
# Kol didžiausia poros |r| >= threshold, pašalina tą poros narį,
# kurio vidutinė |r| su likusiais požymiais didesnė. class nedalyvauja.

remove_correlated_features <- function(data, threshold = 0.95) {

  repeat {
    feature_names <- setdiff(names(data), "class")
    if (length(feature_names) < 2) break

    abs_correlations <- abs(cor(data[feature_names]))
    diag(abs_correlations) <- 0

    if (max(abs_correlations) < threshold) break

    pair <- arrayInd(which.max(abs_correlations), dim(abs_correlations))[1, ]
    pair_names <- feature_names[pair]

    drop <- pair_names[which.max(c(
      mean(abs_correlations[pair_names[1], ]),
      mean(abs_correlations[pair_names[2], ])))]

    cat(sprintf("  šaliname %-18s (|r| = %.3f su %s)\n",
                drop, max(abs_correlations), setdiff(pair_names, drop)))

    data[[drop]] <- NULL
  }

  cat("Likę požymiai:",
      paste(setdiff(names(data), "class"), collapse = ", "), "\n")
  data
}

# UŽDUOTIS 5: z-skaliavimas, z = (x - vidurkis) / standartinis nuokrypis.
# class niekada nėra normalizuojama.

normalize_z_score <- function(data) {
  for (name in setdiff(names(data), "class")) {
    values <- data[[name]]
    data[[name]] <- (values - mean(values)) / sd(values)
  }
  data
}