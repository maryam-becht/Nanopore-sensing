
required_packages <- c(
  "dplyr",
  "purrr",
  "tidyr",
  "tibble",
  "ggplot2",
  "cluster",
  "patchwork",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cluster)
library(patchwork)
library(scales)

set.seed(123)


# ============================================================
# 1. SELECT THE ROOT FOLDER CONTAINING ALL PEPTIDES
# ============================================================

if (interactive() && .Platform$OS.type == "windows") {
  ROOT_FOLDER <- choose.dir(
    default = "C:/Users/mbech/Documents/Peptides",
    caption = "Choisis le dossier qui contient TOUS les dossiers peptides"
  )
} else {
  ROOT_FOLDER <- readline(
    "Entre le chemin complet du dossier contenant tous les peptides : "
  )
}

if (is.na(ROOT_FOLDER) || ROOT_FOLDER == "") {
  stop("Aucun dossier sélectionné.")
}

ROOT_FOLDER <- normalizePath(
  ROOT_FOLDER,
  winslash = "/",
  mustWork = TRUE
)

OUTPUT_FOLDER <- file.path(
  ROOT_FOLDER,
  "Dominant_shape_and_pair_comparison_output"
)

if (!dir.exists(OUTPUT_FOLDER)) {
  dir.create(OUTPUT_FOLDER, recursive = TRUE)
}

message("Root folder: ", ROOT_FOLDER)
message("Output folder: ", OUTPUT_FOLDER)


# ============================================================
# 2. ANALYSIS PARAMETERS
# ============================================================

# Event filters
MIN_DWELL_MS <- 0.5
MAX_DWELL_MS <- 100
MIN_IRES_PCT <- 20
MAX_IRES_PCT <- 85

# Number of points describing each in-event waveform
SHAPE_N_POINTS <- 80

# Global clustering
K_MIN <- 2
K_MAX <- 6

PCA_VARIANCE_TARGET <- 0.95
PCA_MAX_COMPONENTS <- 10

# Balanced training set
MAX_TRAIN_PER_PORE <- 200
MAX_TRAIN_PER_PEPTIDE <- 1000
MAX_TRAIN_TOTAL <- 2500

# Number of real dominant-cluster events displayed per pore
MAX_OVERLAY_PER_PORE <- 10

# Figure settings
PNG_DPI <- 180

# ------------------------------------------------------------
# EDIT THIS TABLE ONLY IF YOUR PEPTIDE NAMES CHANGE
# ------------------------------------------------------------
# These are the complete phospho / non-phospho pairs currently
# present in your peptide dataset.

PHOSPHO_PAIRS <- tibble::tribble(
  ~pair_name,        ~non_phospho,          ~phospho,
  "asyn Y125",       "asynY125",            "asynpY125",
  "tau Thr217",      "tauThr217",           "taupThr217",
  "tau Ser202/Thr205","tauSer202_Thr205",    "AT8"
)


# ============================================================
# 3. FIND ALL EVENT-SHAPE CSV FILES
# ============================================================

shape_files <- list.files(
  path = ROOT_FOLDER,
  pattern = "_nanopore_R_data_event_shapes\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

# Ignore output folders if the script has already been run.
shape_files <- shape_files[
  !grepl(
    "Global_shape_clustering_output|Dominant_shape_and_pair_comparison_output",
    shape_files
  )
]

shape_files <- sort(unique(shape_files))

if (length(shape_files) == 0) {
  stop(
    "Aucun fichier *_nanopore_R_data_event_shapes.csv trouvé sous : ",
    ROOT_FOLDER
  )
}

message("Event-shape CSV files found: ", length(shape_files))


# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

extract_pore_id <- function(file_path) {

  x <- basename(file_path)

  m <- regexpr("P[0-9]+", x, ignore.case = TRUE)

  if (m[1] > 0) {
    return(toupper(regmatches(x, m)))
  }

  tools::file_path_sans_ext(x)
}


extract_peptide_id <- function(file_path, root_folder) {

  parent_dir <- normalizePath(
    dirname(file_path),
    winslash = "/",
    mustWork = FALSE
  )

  root_clean <- normalizePath(
    root_folder,
    winslash = "/",
    mustWork = FALSE
  )

  rel <- sub(
    paste0(
      "^",
      gsub(
        "([.|()\\^{}+$*?]|\\[|\\])",
        "\\\\\\1",
        root_clean
      ),
      "/?"
    ),
    "",
    parent_dir
  )

  rel_parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  rel_parts <- rel_parts[nzchar(rel_parts)]

  if (length(rel_parts) >= 1) {
    rel_parts[1]
  } else {
    basename(parent_dir)
  }
}


safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}


convert_inside_event <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  x_chr <- tolower(trimws(as.character(x)))

  out <- rep(NA, length(x_chr))

  out[x_chr %in% c("1", "true", "t", "yes", "y")] <- TRUE
  out[x_chr %in% c("0", "false", "f", "no", "n")] <- FALSE

  out
}


nearest_centroid_assignment <- function(score_matrix, centers_matrix) {

  n <- nrow(score_matrix)
  k <- nrow(centers_matrix)

  dist_mat <- matrix(
    NA_real_,
    nrow = n,
    ncol = k
  )

  for (j in seq_len(k)) {

    diffs <- sweep(
      score_matrix,
      2,
      centers_matrix[j, ],
      FUN = "-"
    )

    dist_mat[, j] <- rowSums(diffs^2)
  }

  max.col(-dist_mat, ties.method = "first")
}


matrix_rows_to_long <- function(
    mat,
    metadata,
    value_name,
    time_values
) {

  out <- tibble(
    row_index = rep(seq_len(nrow(mat)), each = ncol(mat)),
    normalised_time = rep(time_values, times = nrow(mat)),
    value = as.vector(t(mat))
  )

  out$event_uid <- metadata$event_uid[out$row_index]
  out$peptide <- metadata$peptide[out$row_index]
  out$pore_id <- metadata$pore_id[out$row_index]

  names(out)[names(out) == "value"] <- value_name

  out
}


# ============================================================
# 5. PROCESS ONE EVENT-SHAPE CSV
# ============================================================

process_shape_file <- function(file_path, root_folder) {

  peptide_id <- extract_peptide_id(file_path, root_folder)
  pore_id <- extract_pore_id(file_path)
  source_file <- basename(file_path)

  message("Reading: ", peptide_id, " / ", pore_id)

  dat <- tryCatch(
    utils::read.csv(file_path, check.names = FALSE),
    error = function(e) {
      warning(
        "Could not read ",
        file_path,
        ": ",
        conditionMessage(e)
      )
      return(NULL)
    }
  )

  empty_result <- list(
    metadata = tibble(),
    shape_matrix = matrix(
      numeric(0),
      nrow = 0,
      ncol = SHAPE_N_POINTS
    ),
    display_matrix = matrix(
      numeric(0),
      nrow = 0,
      ncol = SHAPE_N_POINTS
    )
  )

  if (is.null(dat) || nrow(dat) == 0) {
    return(empty_result)
  }

  essential <- c(
    "event_id",
    "time_relative_ms",
    "dwell_total_ms",
    "event_I_over_I0_pct"
  )

  missing_essential <- setdiff(
    essential,
    names(dat)
  )

  if (length(missing_essential) > 0) {

    warning(
      "Skipping ",
      source_file,
      ". Missing columns: ",
      paste(missing_essential, collapse = ", ")
    )

    return(empty_result)
  }

  has_iio_trace <- "I_over_I0_trace_pct" %in% names(dat)

  has_current <- all(
    c("current_pA", "baseline_pA") %in% names(dat)
  )

  if (!has_iio_trace && !has_current) {

    warning(
      "Skipping ",
      source_file,
      ". Need I_over_I0_trace_pct OR current_pA + baseline_pA."
    )

    return(empty_result)
  }

  sh <- tibble(
    event_id = safe_numeric(dat$event_id),
    time_relative_ms = safe_numeric(dat$time_relative_ms),
    dwell_total_ms = safe_numeric(dat$dwell_total_ms),
    event_I_over_I0_pct = safe_numeric(dat$event_I_over_I0_pct)
  )

  if (has_iio_trace) {

    sh$I_over_I0_trace_pct <- safe_numeric(
      dat$I_over_I0_trace_pct
    )

  } else {

    current <- safe_numeric(dat$current_pA)
    baseline <- safe_numeric(dat$baseline_pA)

    sh$I_over_I0_trace_pct <- 100 * current / baseline
  }

  if ("inside_event" %in% names(dat)) {
    sh$inside_event <- convert_inside_event(dat$inside_event)
  } else {
    sh$inside_event <- NA
  }

  # Convert possible 0-1 values to percentages.
  q95_event <- suppressWarnings(
    stats::quantile(
      sh$event_I_over_I0_pct,
      0.95,
      na.rm = TRUE
    )
  )

  if (is.finite(q95_event) && q95_event <= 1.2) {
    sh$event_I_over_I0_pct <- sh$event_I_over_I0_pct * 100
  }

  q95_trace <- suppressWarnings(
    stats::quantile(
      sh$I_over_I0_trace_pct,
      0.95,
      na.rm = TRUE
    )
  )

  if (is.finite(q95_trace) && q95_trace <= 1.2) {
    sh$I_over_I0_trace_pct <- sh$I_over_I0_trace_pct * 100
  }

  sh <- sh |>
    dplyr::filter(
      is.finite(event_id),
      is.finite(time_relative_ms),
      is.finite(dwell_total_ms),
      is.finite(event_I_over_I0_pct),
      is.finite(I_over_I0_trace_pct)
    )

  if (nrow(sh) == 0) {
    return(empty_result)
  }

  event_meta <- sh |>
    dplyr::distinct(
      event_id,
      dwell_total_ms,
      event_I_over_I0_pct
    ) |>
    dplyr::filter(
      dwell_total_ms >= MIN_DWELL_MS,
      dwell_total_ms <= MAX_DWELL_MS,
      event_I_over_I0_pct >= MIN_IRES_PCT,
      event_I_over_I0_pct <= MAX_IRES_PCT
    )

  if (nrow(event_meta) == 0) {
    return(empty_result)
  }

  sh <- sh |>
    dplyr::filter(
      event_id %in% event_meta$event_id
    )

  # Use only the in-event portion.
  if (any(!is.na(sh$inside_event))) {

    sh_inside <- sh |>
      dplyr::filter(inside_event %in% TRUE)

  } else {

    sh_inside <- sh |>
      dplyr::filter(
        time_relative_ms >= 0,
        time_relative_ms <= dwell_total_ms
      )
  }

  if (nrow(sh_inside) == 0) {
    return(empty_result)
  }

  event_splits <- split(
    sh_inside,
    sh_inside$event_id
  )

  grid <- seq(
    0,
    1,
    length.out = SHAPE_N_POINTS
  )

  metadata_list <- list()
  shape_list <- list()
  display_list <- list()

  out_i <- 0L

  for (event_name in names(event_splits)) {

    tr <- event_splits[[event_name]] |>
      dplyr::arrange(time_relative_ms)

    if (nrow(tr) < 6) {
      next
    }

    if (length(unique(tr$time_relative_ms)) < 4) {
      next
    }

    x0 <- min(tr$time_relative_ms, na.rm = TRUE)
    x1 <- max(tr$time_relative_ms, na.rm = TRUE)

    if (
      !is.finite(x0) ||
      !is.finite(x1) ||
      x1 <= x0
    ) {
      next
    }

    x_norm <- (
      tr$time_relative_ms - x0
    ) / (
      x1 - x0
    )

    iio_interp <- tryCatch(
      stats::approx(
        x = x_norm,
        y = tr$I_over_I0_trace_pct,
        xout = grid,
        ties = mean,
        rule = 2
      )$y,
      error = function(e) NULL
    )

    if (
      is.null(iio_interp) ||
      any(!is.finite(iio_interp))
    ) {
      next
    }

    # Shape used for clustering:
    # blockade is positive and amplitude-normalised.
    blockade <- pmax(
      0,
      100 - iio_interp
    )

    amp <- suppressWarnings(
      stats::quantile(
        blockade,
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE
      )
    )

    if (
      !is.finite(amp) ||
      amp <= 1e-8
    ) {
      next
    }

    shape_norm <- blockade / amp

    event_id_num <- safe_numeric(event_name)

    meta_row <- event_meta |>
      dplyr::filter(
        event_id == event_id_num
      ) |>
      dplyr::slice(1)

    if (nrow(meta_row) == 0) {
      next
    }

    out_i <- out_i + 1L

    event_uid <- paste(
      peptide_id,
      pore_id,
      tools::file_path_sans_ext(source_file),
      event_id_num,
      sep = "__"
    )

    metadata_list[[out_i]] <- tibble(
      event_uid = event_uid,
      peptide = peptide_id,
      pore_id = pore_id,
      event_id = event_id_num,
      dwell_total_ms = meta_row$dwell_total_ms[1],
      Ires_pct = meta_row$event_I_over_I0_pct[1],
      source_file = source_file,
      source_path = file_path
    )

    shape_list[[out_i]] <- shape_norm
    display_list[[out_i]] <- iio_interp
  }

  if (out_i == 0L) {
    return(empty_result)
  }

  metadata <- dplyr::bind_rows(
    metadata_list
  )

  shape_matrix <- do.call(
    rbind,
    shape_list
  )

  display_matrix <- do.call(
    rbind,
    display_list
  )

  rownames(shape_matrix) <- metadata$event_uid
  rownames(display_matrix) <- metadata$event_uid

  list(
    metadata = metadata,
    shape_matrix = shape_matrix,
    display_matrix = display_matrix
  )
}


# ============================================================
# 6. READ AND PREPARE ALL FILTERED EVENTS
# ============================================================

processed <- purrr::map(
  shape_files,
  process_shape_file,
  root_folder = ROOT_FOLDER
)

metadata_all <- dplyr::bind_rows(
  purrr::map(
    processed,
    "metadata"
  )
)

valid_processed <- processed[
  vapply(
    processed,
    function(x) nrow(x$metadata) > 0,
    logical(1)
  )
]

if (length(valid_processed) == 0) {
  stop("No valid filtered event shapes found.")
}

shape_matrix_all <- do.call(
  rbind,
  purrr::map(
    valid_processed,
    "shape_matrix"
  )
)

display_matrix_all <- do.call(
  rbind,
  purrr::map(
    valid_processed,
    "display_matrix"
  )
)

if (nrow(metadata_all) < 10) {
  stop(
    "Too few valid filtered events: ",
    nrow(metadata_all)
  )
}

if (
  nrow(shape_matrix_all) != nrow(metadata_all) ||
  nrow(display_matrix_all) != nrow(metadata_all)
) {
  stop(
    "Internal row mismatch between metadata and waveform matrices."
  )
}

message("Filtered events: ", nrow(metadata_all))
message("Peptides: ", length(unique(metadata_all$peptide)))
message(
  "Pore combinations: ",
  nrow(metadata_all |> distinct(peptide, pore_id))
)


# ============================================================
# 7. BUILD BALANCED GLOBAL CLUSTERING TRAINING SET
# ============================================================

training_meta <- metadata_all |>
  dplyr::group_by(peptide, pore_id) |>
  dplyr::group_modify(
    ~ dplyr::slice_sample(
      .x,
      n = min(
        nrow(.x),
        MAX_TRAIN_PER_PORE
      )
    )
  ) |>
  dplyr::ungroup()

training_meta <- training_meta |>
  dplyr::group_by(peptide) |>
  dplyr::group_modify(
    ~ dplyr::slice_sample(
      .x,
      n = min(
        nrow(.x),
        MAX_TRAIN_PER_PEPTIDE
      )
    )
  ) |>
  dplyr::ungroup()

if (nrow(training_meta) > MAX_TRAIN_TOTAL) {
  training_meta <- training_meta |>
    dplyr::slice_sample(
      n = MAX_TRAIN_TOTAL
    )
}

training_idx <- match(
  training_meta$event_uid,
  metadata_all$event_uid
)

training_shape_matrix <- shape_matrix_all[
  training_idx,
  ,
  drop = FALSE
]

message(
  "Balanced training events: ",
  nrow(training_shape_matrix)
)


# ============================================================
# 8. PCA
# ============================================================

pca_fit <- stats::prcomp(
  training_shape_matrix,
  center = TRUE,
  scale. = FALSE
)

var_explained <- (
  pca_fit$sdev^2
) / sum(
  pca_fit$sdev^2
)

cum_var <- cumsum(
  var_explained
)

n_pc <- which(
  cum_var >= PCA_VARIANCE_TARGET
)[1]

if (is.na(n_pc)) {
  n_pc <- length(var_explained)
}

n_pc <- min(
  n_pc,
  PCA_MAX_COMPONENTS,
  ncol(pca_fit$x)
)

n_pc <- max(
  1,
  n_pc
)

training_scores <- pca_fit$x[
  ,
  seq_len(n_pc),
  drop = FALSE
]

message(
  "PCA components retained: ",
  n_pc,
  " (",
  sprintf(
    "%.1f%% cumulative variance",
    100 * cum_var[n_pc]
  ),
  ")"
)


# ============================================================
# 9. CHOOSE GLOBAL NUMBER OF CLUSTERS
# ============================================================

max_valid_k <- min(
  K_MAX,
  nrow(training_scores) - 1
)

k_candidates <- seq.int(
  K_MIN,
  max_valid_k
)

if (length(k_candidates) == 0) {
  stop("Not enough training events to test clustering.")
}

training_dist <- stats::dist(
  training_scores
)

k_results <- purrr::map_dfr(
  k_candidates,
  function(k) {

    km_try <- tryCatch(
      stats::kmeans(
        training_scores,
        centers = k,
        nstart = 30,
        iter.max = 200
      ),
      error = function(e) NULL
    )

    if (is.null(km_try)) {
      return(
        tibble(
          k = k,
          average_silhouette = NA_real_
        )
      )
    }

    sil <- cluster::silhouette(
      km_try$cluster,
      training_dist
    )

    tibble(
      k = k,
      average_silhouette = mean(
        sil[, "sil_width"],
        na.rm = TRUE
      )
    )
  }
)

valid_k_results <- k_results |>
  dplyr::filter(
    is.finite(average_silhouette)
  )

if (nrow(valid_k_results) == 0) {
  stop(
    "Could not obtain a valid silhouette score."
  )
}

best_k <- valid_k_results |>
  dplyr::arrange(
    dplyr::desc(average_silhouette),
    k
  ) |>
  dplyr::slice(1) |>
  dplyr::pull(k)

message(
  "Selected global number of clusters: ",
  best_k
)

write.csv(
  k_results,
  file.path(
    OUTPUT_FOLDER,
    "00_cluster_number_silhouette_scores.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. FIT FINAL GLOBAL CLUSTER MODEL
# ============================================================

final_kmeans <- stats::kmeans(
  training_scores,
  centers = best_k,
  nstart = 100,
  iter.max = 300
)

all_scores_full <- stats::predict(
  pca_fit,
  newdata = shape_matrix_all
)

all_scores <- all_scores_full[
  ,
  seq_len(n_pc),
  drop = FALSE
]

raw_cluster <- nearest_centroid_assignment(
  score_matrix = all_scores,
  centers_matrix = final_kmeans$centers
)

metadata_all$cluster_raw <- raw_cluster

# Cluster 1 = globally most abundant cluster.
cluster_order <- metadata_all |>
  dplyr::count(
    cluster_raw,
    name = "N"
  ) |>
  dplyr::arrange(
    dplyr::desc(N),
    cluster_raw
  ) |>
  dplyr::mutate(
    cluster = dplyr::row_number()
  )

relabel_map <- setNames(
  cluster_order$cluster,
  cluster_order$cluster_raw
)

metadata_all$cluster <- as.integer(
  relabel_map[
    as.character(
      metadata_all$cluster_raw
    )
  ]
)

metadata_all$cluster_label <- factor(
  paste0(
    "Cluster ",
    metadata_all$cluster
  ),
  levels = paste0(
    "Cluster ",
    seq_len(best_k)
  )
)

all_cluster_levels <- levels(
  metadata_all$cluster_label
)


# ============================================================
# 11. SAVE EVENT ASSIGNMENTS
# ============================================================

write.csv(
  metadata_all |>
    dplyr::select(
      peptide,
      pore_id,
      event_id,
      dwell_total_ms,
      Ires_pct,
      cluster,
      cluster_label,
      source_file,
      source_path,
      event_uid
    ),
  file.path(
    OUTPUT_FOLDER,
    "01_all_filtered_events_with_global_clusters.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. CORRECT CLUSTER PROPORTIONS PER PORE
# ============================================================

# IMPORTANT:
# Build cluster rows only from REAL peptide/pore combinations.
# This avoids creating phantom pores.

actual_pores <- metadata_all |>
  dplyr::distinct(
    peptide,
    pore_id
  )

pore_cluster_counts <- metadata_all |>
  dplyr::count(
    peptide,
    pore_id,
    cluster_label,
    name = "N_events"
  )

pore_cluster_grid <- tidyr::crossing(
  actual_pores,
  cluster_label = factor(
    all_cluster_levels,
    levels = all_cluster_levels
  )
)

pore_props <- pore_cluster_grid |>
  dplyr::left_join(
    pore_cluster_counts,
    by = c(
      "peptide",
      "pore_id",
      "cluster_label"
    )
  ) |>
  dplyr::mutate(
    N_events = dplyr::coalesce(
      N_events,
      0L
    )
  ) |>
  dplyr::group_by(
    peptide,
    pore_id
  ) |>
  dplyr::mutate(
    N_total_filtered = sum(N_events),
    proportion = N_events / N_total_filtered
  ) |>
  dplyr::ungroup()

write.csv(
  pore_props,
  file.path(
    OUTPUT_FOLDER,
    "02_cluster_proportions_by_pore.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. DOMINANT SHAPE CLUSTER PER PORE
# ============================================================

pore_dominant <- pore_props |>
  dplyr::group_by(
    peptide,
    pore_id
  ) |>
  dplyr::arrange(
    dplyr::desc(proportion),
    cluster_label,
    .by_group = TRUE
  ) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::rename(
    dominant_cluster_label = cluster_label,
    dominant_N_events = N_events,
    dominant_proportion = proportion
  )

pore_dominant <- pore_dominant |>
  dplyr::mutate(
    dominant_cluster = as.integer(
      sub(
        "Cluster ",
        "",
        as.character(
          dominant_cluster_label
        )
      )
    )
  )

write.csv(
  pore_dominant,
  file.path(
    OUTPUT_FOLDER,
    "03_dominant_cluster_by_pore.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. BUILD DISPLAY TABLES
# ============================================================

norm_time <- seq(
  0,
  1,
  length.out = SHAPE_N_POINTS
)

display_long <- matrix_rows_to_long(
  display_matrix_all,
  metadata_all,
  "I_over_I0_pct",
  norm_time
)

shape_long <- matrix_rows_to_long(
  shape_matrix_all,
  metadata_all,
  "normalised_blockade",
  norm_time
)

display_long <- display_long |>
  dplyr::left_join(
    metadata_all |>
      dplyr::select(
        event_uid,
        cluster,
        cluster_label,
        dwell_total_ms,
        Ires_pct
      ),
    by = "event_uid"
  )

shape_long <- shape_long |>
  dplyr::left_join(
    metadata_all |>
      dplyr::select(
        event_uid,
        cluster,
        cluster_label,
        dwell_total_ms,
        Ires_pct
      ),
    by = "event_uid"
  )


# ============================================================
# 15. PORE DOMINANT MEDIAN SHAPES
# ============================================================

pore_dominant_events <- metadata_all |>
  dplyr::inner_join(
    pore_dominant |>
      dplyr::select(
        peptide,
        pore_id,
        dominant_cluster
      ),
    by = c(
      "peptide",
      "pore_id"
    )
  ) |>
  dplyr::filter(
    cluster == dominant_cluster
  )

pore_median_display <- display_long |>
  dplyr::inner_join(
    pore_dominant_events |>
      dplyr::select(event_uid),
    by = "event_uid"
  ) |>
  dplyr::group_by(
    peptide,
    pore_id,
    normalised_time
  ) |>
  dplyr::summarise(
    median_I_over_I0_pct = stats::median(
      I_over_I0_pct,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

pore_median_shape <- shape_long |>
  dplyr::inner_join(
    pore_dominant_events |>
      dplyr::select(event_uid),
    by = "event_uid"
  ) |>
  dplyr::group_by(
    peptide,
    pore_id,
    normalised_time
  ) |>
  dplyr::summarise(
    median_normalised_blockade = stats::median(
      normalised_blockade,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ============================================================
# 16. SELECT REAL DOMINANT EVENTS FOR PORE OVERLAYS
# ============================================================

overlay_event_uids <- character(0)

for (i in seq_len(nrow(pore_dominant))) {

  pep <- pore_dominant$peptide[i]
  por <- pore_dominant$pore_id[i]
  clu <- pore_dominant$dominant_cluster[i]

  idx <- which(
    metadata_all$peptide == pep &
    metadata_all$pore_id == por &
    metadata_all$cluster == clu
  )

  if (length(idx) == 0) {
    next
  }

  mat <- shape_matrix_all[
    idx,
    ,
    drop = FALSE
  ]

  med <- apply(
    mat,
    2,
    stats::median,
    na.rm = TRUE
  )

  d <- apply(
    mat,
    1,
    function(z) {
      sqrt(
        mean(
          (z - med)^2,
          na.rm = TRUE
        )
      )
    }
  )

  keep_n <- min(
    MAX_OVERLAY_PER_PORE,
    length(idx)
  )

  keep_local <- order(d)[
    seq_len(keep_n)
  ]

  overlay_event_uids <- c(
    overlay_event_uids,
    metadata_all$event_uid[
      idx[keep_local]
    ]
  )
}

overlay_event_uids <- unique(
  overlay_event_uids
)

pore_overlay_long <- display_long |>
  dplyr::filter(
    event_uid %in% overlay_event_uids
  )


# ============================================================
# 17. FIGURE: MOST FREQUENT SHAPE IN EACH PORE
# ============================================================

pore_panel_labels <- pore_dominant |>
  dplyr::mutate(
    panel = paste0(
      peptide,
      " | ",
      pore_id,
      "\n",
      dominant_cluster_label,
      " = ",
      scales::percent(
        dominant_proportion,
        accuracy = 1
      ),
      " (n=",
      dominant_N_events,
      ")"
    )
  )

pore_overlay_plot_data <- pore_overlay_long |>
  dplyr::left_join(
    pore_panel_labels |>
      dplyr::select(
        peptide,
        pore_id,
        panel
      ),
    by = c(
      "peptide",
      "pore_id"
    )
  )

pore_median_plot_data <- pore_median_display |>
  dplyr::left_join(
    pore_panel_labels |>
      dplyr::select(
        peptide,
        pore_id,
        panel
      ),
    by = c(
      "peptide",
      "pore_id"
    )
  )

p_pore_shapes <- ggplot2::ggplot() +

  ggplot2::geom_line(
    data = pore_overlay_plot_data,
    ggplot2::aes(
      x = normalised_time,
      y = I_over_I0_pct,
      group = event_uid
    ),
    colour = "grey78",
    linewidth = 0.25,
    alpha = 0.65
  ) +

  ggplot2::geom_line(
    data = pore_median_plot_data,
    ggplot2::aes(
      x = normalised_time,
      y = median_I_over_I0_pct
    ),
    colour = "black",
    linewidth = 0.80
  ) +

  ggplot2::facet_wrap(
    ~ panel,
    ncol = 4,
    scales = "free_y"
  ) +

  ggplot2::labs(
    x = "Normalised event time",
    y = expression(I/I[0] ~ " (%)"),
    title = "Most frequent event-shape population in each pore",
    subtitle = "Grey = real events from the dominant global cluster; black = pore median waveform"
  ) +

  ggplot2::theme_classic(
    base_size = 8
  ) +

  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 7
    ),
    legend.position = "none"
  )

ggplot2::ggsave(
  file.path(
    OUTPUT_FOLDER,
    "04_most_frequent_shape_in_each_pore.png"
  ),
  p_pore_shapes,
  width = 13.0,
  height = max(
    7.0,
    2.4 * ceiling(
      nrow(pore_dominant) / 4
    )
  ),
  dpi = PNG_DPI,
  bg = "white",
  limitsize = FALSE
)


# ============================================================
# 18. PORE-BALANCED CLUSTER PROPORTIONS PER PEPTIDE
# ============================================================

peptide_props <- pore_props |>
  dplyr::group_by(
    peptide,
    cluster_label
  ) |>
  dplyr::summarise(
    mean_proportion = mean(
      proportion,
      na.rm = TRUE
    ),
    sd_proportion = stats::sd(
      proportion,
      na.rm = TRUE
    ),
    N_pores = dplyr::n_distinct(
      pore_id
    ),
    .groups = "drop"
  )

write.csv(
  peptide_props,
  file.path(
    OUTPUT_FOLDER,
    "05_cluster_proportions_by_peptide_pore_balanced.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. DOMINANT CLUSTER PER PEPTIDE
# ============================================================

peptide_dominant <- peptide_props |>
  dplyr::group_by(
    peptide
  ) |>
  dplyr::arrange(
    dplyr::desc(mean_proportion),
    cluster_label,
    .by_group = TRUE
  ) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::rename(
    dominant_cluster_label = cluster_label,
    dominant_mean_proportion = mean_proportion,
    dominant_sd_proportion = sd_proportion
  ) |>
  dplyr::mutate(
    dominant_cluster = as.integer(
      sub(
        "Cluster ",
        "",
        as.character(
          dominant_cluster_label
        )
      )
    )
  )

# How many pores have the same dominant cluster as the peptide?
pore_cluster_agreement <- pore_dominant |>
  dplyr::left_join(
    peptide_dominant |>
      dplyr::select(
        peptide,
        peptide_dominant_cluster = dominant_cluster
      ),
    by = "peptide"
  ) |>
  dplyr::group_by(
    peptide
  ) |>
  dplyr::summarise(
    N_pores_total = dplyr::n(),
    N_pores_matching_peptide_dominant = sum(
      dominant_cluster == peptide_dominant_cluster
    ),
    .groups = "drop"
  )

peptide_dominant <- peptide_dominant |>
  dplyr::left_join(
    pore_cluster_agreement,
    by = "peptide"
  )

write.csv(
  peptide_dominant,
  file.path(
    OUTPUT_FOLDER,
    "06_dominant_cluster_by_peptide.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. BUILD PORE-BALANCED PEPTIDE CONSENSUS SHAPES
# ============================================================

# For each peptide, use events belonging to the PEPTIDE'S dominant
# global cluster. First calculate the median within each pore,
# then take the median across pores so every pore contributes equally.

peptide_dom_event_meta <- metadata_all |>
  dplyr::left_join(
    peptide_dominant |>
      dplyr::select(
        peptide,
        peptide_dominant_cluster = dominant_cluster
      ),
    by = "peptide"
  ) |>
  dplyr::filter(
    cluster == peptide_dominant_cluster
  )

# Per-pore median in actual I/I0 units
peptide_dom_pore_display <- display_long |>
  dplyr::inner_join(
    peptide_dom_event_meta |>
      dplyr::select(
        event_uid,
        peptide,
        pore_id
      ),
    by = c(
      "event_uid",
      "peptide",
      "pore_id"
    )
  ) |>
  dplyr::group_by(
    peptide,
    pore_id,
    normalised_time
  ) |>
  dplyr::summarise(
    pore_median_I_over_I0_pct = stats::median(
      I_over_I0_pct,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Peptide consensus = median of pore medians
peptide_consensus_display <- peptide_dom_pore_display |>
  dplyr::group_by(
    peptide,
    normalised_time
  ) |>
  dplyr::summarise(
    consensus_I_over_I0_pct = stats::median(
      pore_median_I_over_I0_pct,
      na.rm = TRUE
    ),
    Q25_pore = stats::quantile(
      pore_median_I_over_I0_pct,
      0.25,
      na.rm = TRUE
    ),
    Q75_pore = stats::quantile(
      pore_median_I_over_I0_pct,
      0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Same calculation for the amplitude-normalised waveform shape.
peptide_dom_pore_shape <- shape_long |>
  dplyr::inner_join(
    peptide_dom_event_meta |>
      dplyr::select(
        event_uid,
        peptide,
        pore_id
      ),
    by = c(
      "event_uid",
      "peptide",
      "pore_id"
    )
  ) |>
  dplyr::group_by(
    peptide,
    pore_id,
    normalised_time
  ) |>
  dplyr::summarise(
    pore_median_normalised_blockade = stats::median(
      normalised_blockade,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

peptide_consensus_shape <- peptide_dom_pore_shape |>
  dplyr::group_by(
    peptide,
    normalised_time
  ) |>
  dplyr::summarise(
    consensus_normalised_blockade = stats::median(
      pore_median_normalised_blockade,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  peptide_consensus_display,
  file.path(
    OUTPUT_FOLDER,
    "07_peptide_dominant_consensus_waveforms_actual_Ires.csv"
  ),
  row.names = FALSE
)

write.csv(
  peptide_consensus_shape,
  file.path(
    OUTPUT_FOLDER,
    "08_peptide_dominant_consensus_waveforms_normalised_shape.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. FIGURE: MOST FREQUENT SHAPE IN EACH PEPTIDE
# ============================================================

peptide_labels <- peptide_dominant |>
  dplyr::mutate(
    panel = paste0(
      peptide,
      "\n",
      dominant_cluster_label,
      " = ",
      scales::percent(
        dominant_mean_proportion,
        accuracy = 1
      ),
      " | ",
      N_pores_matching_peptide_dominant,
      "/",
      N_pores_total,
      " pores agree"
    )
  )

peptide_pore_plot_data <- peptide_dom_pore_display |>
  dplyr::left_join(
    peptide_labels |>
      dplyr::select(
        peptide,
        panel
      ),
    by = "peptide"
  )

peptide_consensus_plot_data <- peptide_consensus_display |>
  dplyr::left_join(
    peptide_labels |>
      dplyr::select(
        peptide,
        panel
      ),
    by = "peptide"
  )

p_peptide_shapes <- ggplot2::ggplot() +

  ggplot2::geom_line(
    data = peptide_pore_plot_data,
    ggplot2::aes(
      x = normalised_time,
      y = pore_median_I_over_I0_pct,
      group = pore_id
    ),
    colour = "grey72",
    linewidth = 0.45,
    alpha = 0.80
  ) +

  ggplot2::geom_line(
    data = peptide_consensus_plot_data,
    ggplot2::aes(
      x = normalised_time,
      y = consensus_I_over_I0_pct
    ),
    colour = "black",
    linewidth = 0.95
  ) +

  ggplot2::facet_wrap(
    ~ panel,
    ncol = 3,
    scales = "free_y"
  ) +

  ggplot2::labs(
    x = "Normalised event time",
    y = expression(I/I[0] ~ " (%)"),
    title = "Most frequent event shape for each peptide",
    subtitle = "Grey = pore-level median shapes; black = pore-balanced peptide consensus"
  ) +

  ggplot2::theme_classic(
    base_size = 9
  ) +

  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 8
    ),
    legend.position = "none"
  )

ggplot2::ggsave(
  file.path(
    OUTPUT_FOLDER,
    "09_most_frequent_shape_in_each_peptide.png"
  ),
  p_peptide_shapes,
  width = 11.5,
  height = max(
    6.0,
    3.0 * ceiling(
      nrow(peptide_dominant) / 3
    )
  ),
  dpi = PNG_DPI,
  bg = "white",
  limitsize = FALSE
)


# ============================================================
# 22. FIGURE: ALL PEPTIDE CONSENSUS SHAPES OVERLAID
# ============================================================

p_all_peptides <- ggplot2::ggplot(
  peptide_consensus_shape,
  ggplot2::aes(
    x = normalised_time,
    y = consensus_normalised_blockade,
    colour = peptide,
    group = peptide
  )
) +

  ggplot2::geom_line(
    linewidth = 0.85
  ) +

  ggplot2::labs(
    x = "Normalised event time",
    y = "Normalised blockade shape",
    colour = "Peptide",
    title = "Comparison of dominant event shapes across peptides",
    subtitle = "Amplitude-normalised shapes: comparison focuses on waveform form"
  ) +

  ggplot2::theme_classic(
    base_size = 9
  ) +

  ggplot2::theme(
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(
    OUTPUT_FOLDER,
    "10_all_peptide_dominant_shapes_overlay_normalised.png"
  ),
  p_all_peptides,
  width = 9.5,
  height = 5.8,
  dpi = PNG_DPI,
  bg = "white"
)


# ============================================================
# 23. PAIR DEFINITIONS PRESENT IN THE CURRENT DATA
# ============================================================

available_peptides <- unique(
  metadata_all$peptide
)

pairs_present <- PHOSPHO_PAIRS |>
  dplyr::filter(
    non_phospho %in% available_peptides,
    phospho %in% available_peptides
  )

if (nrow(pairs_present) == 0) {

  warning(
    "No complete phospho/non-phospho pair from PHOSPHO_PAIRS was found."
  )

} else {

  # ==========================================================
  # 24. BUILD PAIR WAVEFORM DATA
  # ==========================================================

  pair_shape_display <- pairs_present |>
    tidyr::pivot_longer(
      cols = c(
        non_phospho,
        phospho
      ),
      names_to = "status",
      values_to = "peptide"
    ) |>
    dplyr::mutate(
      status = dplyr::recode(
        status,
        non_phospho = "Non-phosphorylated",
        phospho = "Phosphorylated"
      )
    ) |>
    dplyr::left_join(
      peptide_consensus_display,
      by = "peptide"
    )

  pair_shape_normalised <- pairs_present |>
    tidyr::pivot_longer(
      cols = c(
        non_phospho,
        phospho
      ),
      names_to = "status",
      values_to = "peptide"
    ) |>
    dplyr::mutate(
      status = dplyr::recode(
        status,
        non_phospho = "Non-phosphorylated",
        phospho = "Phosphorylated"
      )
    ) |>
    dplyr::left_join(
      peptide_consensus_shape,
      by = "peptide"
    )


  # ==========================================================
  # 25. PHOSPHO PAIR SHAPE DISTANCES
  # ==========================================================

  pair_distance_list <- vector(
    "list",
    nrow(pairs_present)
  )

  for (i in seq_len(nrow(pairs_present))) {

    pair_i <- pairs_present[i, ]

    np_name <- pair_i$non_phospho
    p_name <- pair_i$phospho

    np_shape <- peptide_consensus_shape |>
      dplyr::filter(
        peptide == np_name
      ) |>
      dplyr::arrange(
        normalised_time
      )

    p_shape <- peptide_consensus_shape |>
      dplyr::filter(
        peptide == p_name
      ) |>
      dplyr::arrange(
        normalised_time
      )

    np_actual <- peptide_consensus_display |>
      dplyr::filter(
        peptide == np_name
      ) |>
      dplyr::arrange(
        normalised_time
      )

    p_actual <- peptide_consensus_display |>
      dplyr::filter(
        peptide == p_name
      ) |>
      dplyr::arrange(
        normalised_time
      )

    if (
      nrow(np_shape) == SHAPE_N_POINTS &&
      nrow(p_shape) == SHAPE_N_POINTS
    ) {

      shape_rmse <- sqrt(
        mean(
          (
            np_shape$consensus_normalised_blockade -
              p_shape$consensus_normalised_blockade
          )^2,
          na.rm = TRUE
        )
      )

    } else {
      shape_rmse <- NA_real_
    }

    if (
      nrow(np_actual) == SHAPE_N_POINTS &&
      nrow(p_actual) == SHAPE_N_POINTS
    ) {

      actual_rmse <- sqrt(
        mean(
          (
            np_actual$consensus_I_over_I0_pct -
              p_actual$consensus_I_over_I0_pct
          )^2,
          na.rm = TRUE
        )
      )

    } else {
      actual_rmse <- NA_real_
    }

    pair_distance_list[[i]] <- tibble(
      pair_name = pair_i$pair_name,
      non_phospho = np_name,
      phospho = p_name,
      normalised_shape_RMSE = shape_rmse,
      actual_Ires_profile_RMSE_pct = actual_rmse
    )
  }

  pair_distances <- dplyr::bind_rows(
    pair_distance_list
  )

  write.csv(
    pair_distances,
    file.path(
      OUTPUT_FOLDER,
      "11_phospho_pair_shape_distances.csv"
    ),
    row.names = FALSE
  )


  # ==========================================================
  # 26. FIGURE: PHOSPHO/NON-PHOSPHO WAVEFORM COMPARISON
  # ==========================================================

  p_pair_waveforms <- ggplot2::ggplot(
    pair_shape_display,
    ggplot2::aes(
      x = normalised_time,
      y = consensus_I_over_I0_pct,
      colour = status,
      linetype = status,
      group = status
    )
  ) +

    ggplot2::geom_line(
      linewidth = 0.95
    ) +

    ggplot2::facet_wrap(
      ~ pair_name,
      ncol = 1,
      scales = "free_y"
    ) +

    ggplot2::labs(
      x = "Normalised event time",
      y = expression(I/I[0] ~ " (%)"),
      colour = NULL,
      linetype = NULL,
      title = "Dominant event-shape comparison: phospho vs non-phospho"
    ) +

    ggplot2::theme_classic(
      base_size = 9
    ) +

    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      )
    )


  # ==========================================================
  # 27. PAIR CLUSTER PROPORTIONS
  # ==========================================================

  pair_peptide_info <- pairs_present |>
    tidyr::pivot_longer(
      cols = c(
        non_phospho,
        phospho
      ),
      names_to = "status",
      values_to = "peptide"
    ) |>
    dplyr::mutate(
      status = dplyr::recode(
        status,
        non_phospho = "Non-phosphorylated",
        phospho = "Phosphorylated"
      )
    )

  pair_props <- pair_peptide_info |>
    dplyr::left_join(
      peptide_props,
      by = "peptide"
    )

  p_pair_props <- ggplot2::ggplot(
    pair_props,
    ggplot2::aes(
      x = status,
      y = mean_proportion,
      fill = cluster_label
    )
  ) +

    ggplot2::geom_col(
      width = 0.68,
      colour = "white",
      linewidth = 0.20
    ) +

    ggplot2::facet_wrap(
      ~ pair_name,
      ncol = 1
    ) +

    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(
        accuracy = 1
      ),
      expand = ggplot2::expansion(
        mult = c(0, 0.02)
      )
    ) +

    ggplot2::labs(
      x = NULL,
      y = "Pore-balanced cluster proportion",
      fill = "Shape cluster",
      title = "Distribution of waveform-shape clusters"
    ) +

    ggplot2::theme_classic(
      base_size = 9
    ) +

    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 20,
        hjust = 1
      )
    )


  # ==========================================================
  # 28. COMBINED PAIR FIGURE:
  #     SHAPES + CLUSTER PROPORTIONS TOGETHER
  # ==========================================================

  p_pair_combined <- (
    p_pair_waveforms +
      p_pair_props +
      patchwork::plot_layout(
        widths = c(1.65, 1)
      )
  ) +
    patchwork::plot_annotation(
      title = "Phosphorylation comparison: dominant waveform shape and event-shape populations"
    )

  ggplot2::ggsave(
    file.path(
      OUTPUT_FOLDER,
      "12_phospho_pairs_shapes_AND_cluster_proportions.png"
    ),
    p_pair_combined,
    width = 12.5,
    height = max(
      7.2,
      2.6 * nrow(pairs_present) + 1.5
    ),
    dpi = PNG_DPI,
    bg = "white",
    limitsize = FALSE
  )


  # ==========================================================
  # 29. NORMALISED PAIR SHAPE FIGURE
  # ==========================================================

  p_pair_normalised <- ggplot2::ggplot(
    pair_shape_normalised,
    ggplot2::aes(
      x = normalised_time,
      y = consensus_normalised_blockade,
      colour = status,
      linetype = status,
      group = status
    )
  ) +

    ggplot2::geom_line(
      linewidth = 0.95
    ) +

    ggplot2::facet_wrap(
      ~ pair_name,
      ncol = 1
    ) +

    ggplot2::labs(
      x = "Normalised event time",
      y = "Normalised blockade shape",
      colour = NULL,
      linetype = NULL,
      title = "Phospho vs non-phospho: waveform FORM only",
      subtitle = "Amplitude normalised to focus on differences in shape"
    ) +

    ggplot2::theme_classic(
      base_size = 9
    ) +

    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      )
    )

  ggplot2::ggsave(
    file.path(
      OUTPUT_FOLDER,
      "13_phospho_pairs_normalised_shape_only.png"
    ),
    p_pair_normalised,
    width = 8.5,
    height = max(
      6.5,
      2.4 * nrow(pairs_present) + 1.2
    ),
    dpi = PNG_DPI,
    bg = "white",
    limitsize = FALSE
  )
}


# ============================================================
# 30. ALL-PEPTIDE PAIRWISE DOMINANT-SHAPE DISTANCES
# ============================================================

peptide_names <- sort(
  unique(
    peptide_consensus_shape$peptide
  )
)

shape_wide <- peptide_consensus_shape |>
  dplyr::select(
    peptide,
    normalised_time,
    consensus_normalised_blockade
  ) |>
  tidyr::pivot_wider(
    names_from = peptide,
    values_from = consensus_normalised_blockade
  ) |>
  dplyr::arrange(
    normalised_time
  )

distance_rows <- list()
counter <- 0L

for (i in seq_along(peptide_names)) {

  for (j in seq_along(peptide_names)) {

    pep1 <- peptide_names[i]
    pep2 <- peptide_names[j]

    v1 <- shape_wide[[pep1]]
    v2 <- shape_wide[[pep2]]

    rmse <- sqrt(
      mean(
        (v1 - v2)^2,
        na.rm = TRUE
      )
    )

    counter <- counter + 1L

    distance_rows[[counter]] <- tibble(
      peptide_1 = pep1,
      peptide_2 = pep2,
      shape_RMSE = rmse
    )
  }
}

all_pairwise_distances <- dplyr::bind_rows(
  distance_rows
)

write.csv(
  all_pairwise_distances,
  file.path(
    OUTPUT_FOLDER,
    "14_all_peptide_pairwise_dominant_shape_distances.csv"
  ),
  row.names = FALSE
)

p_distance_heatmap <- ggplot2::ggplot(
  all_pairwise_distances,
  ggplot2::aes(
    x = peptide_1,
    y = peptide_2,
    fill = shape_RMSE
  )
) +

  ggplot2::geom_tile() +

  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.3f",
        shape_RMSE
      )
    ),
    size = 2.7
  ) +

  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "Shape RMSE",
    title = "Difference between dominant waveform shapes across peptides",
    subtitle = "Lower values = more similar dominant waveform form"
  ) +

  ggplot2::theme_classic(
    base_size = 8
  ) +

  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

ggplot2::ggsave(
  file.path(
    OUTPUT_FOLDER,
    "15_all_peptide_dominant_shape_distance_heatmap.png"
  ),
  p_distance_heatmap,
  width = 8.5,
  height = 7.2,
  dpi = PNG_DPI,
  bg = "white"
)


# ============================================================
# 31. README / SUMMARY
# ============================================================

readme_lines <- c(
  paste0("Root folder: ", ROOT_FOLDER),
  paste0("Event-shape CSV files: ", length(shape_files)),
  paste0("Filtered events: ", nrow(metadata_all)),
  paste0(
    "Filter: dwell >= ",
    MIN_DWELL_MS,
    " ms and Ires = ",
    MIN_IRES_PCT,
    "-",
    MAX_IRES_PCT,
    "%"
  ),
  paste0("Global clusters selected: ", best_k),
  "",
  "HOW TO READ THE OUTPUTS",
  "",
  "04_most_frequent_shape_in_each_pore.png",
  "  - Each panel is one pore.",
  "  - Grey lines are real events belonging to that pore's dominant GLOBAL shape cluster.",
  "  - Black line is the median waveform for that dominant pore population.",
  "",
  "09_most_frequent_shape_in_each_peptide.png",
  "  - Each panel is one peptide.",
  "  - Grey lines are pore-level median waveforms.",
  "  - Black line is the pore-balanced peptide consensus waveform.",
  "",
  "10_all_peptide_dominant_shapes_overlay_normalised.png",
  "  - Direct comparison of the dominant waveform FORM across all peptides.",
  "",
  "12_phospho_pairs_shapes_AND_cluster_proportions.png",
  "  - LEFT: dominant peptide waveform for each phospho/non-phospho pair.",
  "  - RIGHT: pore-balanced proportion of each global shape cluster.",
  "  - This is the main figure for asking whether phosphorylation changes",
  "    both the characteristic waveform and the distribution of event-shape populations.",
  "",
  "13_phospho_pairs_normalised_shape_only.png",
  "  - Compares waveform form after amplitude normalisation.",
  "",
  "15_all_peptide_dominant_shape_distance_heatmap.png",
  "  - Quantifies similarity between peptide dominant shapes.",
  "  - Lower RMSE = more similar shape.",
  "",
  "IMPORTANT",
  "  Global clusters are shared across all peptides.",
  "  Dominant shape per pore/peptide means the most frequent GLOBAL cluster",
  "  in that pore/peptide, followed by a median waveform calculated within that population.",
  "  Peptide-level proportions and consensus shapes are pore-balanced."
)

writeLines(
  readme_lines,
  file.path(
    OUTPUT_FOLDER,
    "README_dominant_shape_and_pair_comparison.txt"
  )
)

message("")
message("DONE.")
message("Output folder: ", OUTPUT_FOLDER)
message("")
message("Main figures to inspect first:")
message("  04_most_frequent_shape_in_each_pore.png")
message("  09_most_frequent_shape_in_each_peptide.png")
message("  12_phospho_pairs_shapes_AND_cluster_proportions.png")
message("  15_all_peptide_dominant_shape_distance_heatmap.png")
