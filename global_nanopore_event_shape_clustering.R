# ============================================================
# GLOBAL NANOPORE EVENT-SHAPE CLUSTERING ACROSS ALL PEPTIDES
# ============================================================
#
# PURPOSE
#   1) Read every *_nanopore_R_data_event_shapes.csv under a root folder.
#   2) Keep events with:
#        dwell time >= 0.5 ms
#        Ires between 20 and 85 %
#   3) Compare waveform SHAPE only:
#        - use the in-event trace
#        - normalise event time to 0-1
#        - normalise blockade amplitude for clustering only
#   4) Build GLOBAL clusters shared by all peptides/pores.
#   5) Choose the number of clusters objectively using silhouette score.
#   6) Assign every filtered event to one global cluster.
#   7) Generate:
#        - representative shapes for each global cluster
#        - Ires vs dwell scatter coloured by cluster
#        - cluster proportions for each peptide
#        - CSV tables with assignments and summaries
#
# IMPORTANT
#   Clustering is performed on ALL filtered event waveforms, not only on
#   pre-selected representative events. Representative shapes/events are
#   extracted AFTER clustering.
#
# EXPECTED FOLDER STRUCTURE
#   Peptides/
#     tauThr217/
#       tauThr217_P2_nanopore_R_data_event_shapes.csv
#       tauThr217_P4_nanopore_R_data_event_shapes.csv
#     taupThr217/
#       taupThr217_P1_nanopore_R_data_event_shapes.csv
#     AT8/
#       AT8_P1_nanopore_R_data_event_shapes.csv
#     ...
#
# OUTPUT FOLDER
#   Peptides/Global_shape_clustering_output/
#
# ============================================================

# ============================================================
# 0. PACKAGES
# ============================================================

required_packages <- c(
  "dplyr",
  "purrr",
  "tidyr",
  "tibble",
  "ggplot2",
  "cluster",
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

OUTPUT_FOLDER <- file.path(ROOT_FOLDER, "Global_shape_clustering_output")
if (!dir.exists(OUTPUT_FOLDER)) {
  dir.create(OUTPUT_FOLDER, recursive = TRUE)
}

message("Root folder: ", ROOT_FOLDER)
message("Output folder: ", OUTPUT_FOLDER)


# ============================================================
# 2. ANALYSIS PARAMETERS
# ============================================================

# Filters announced for the representative-event analysis
MIN_DWELL_MS <- 0.5
MAX_DWELL_MS <- 100
MIN_IRES_PCT <- 20
MAX_IRES_PCT <- 85

# Number of points used to describe each event shape
SHAPE_N_POINTS <- 80

# K values tested for global clustering
K_MIN <- 2
K_MAX <- 6

# PCA variance retained before clustering
PCA_VARIANCE_TARGET <- 0.95
PCA_MAX_COMPONENTS <- 10

# Balanced training set: prevents one very event-rich pore from
# defining the global clusters by itself.
MAX_TRAIN_PER_PORE <- 200
MAX_TRAIN_PER_PEPTIDE <- 1000
MAX_TRAIN_TOTAL <- 2500

# Number of real event traces overlaid in each cluster-shape panel
MAX_OVERLAY_PER_CLUSTER <- 20

# Plot settings
PNG_DPI <- 180


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

# Ignore anything accidentally placed inside the global output folder.
shape_files <- shape_files[
  !grepl("Global_shape_clustering_output", shape_files, fixed = TRUE)
]

shape_files <- sort(unique(shape_files))

if (length(shape_files) == 0) {
  stop(
    "Aucun fichier *_nanopore_R_data_event_shapes.csv trouvé sous : ",
    ROOT_FOLDER
  )
}

message("Event-shape CSV files found: ", length(shape_files))
print(shape_files)


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
    paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", root_clean), "/?"),
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

  dist_mat <- matrix(NA_real_, nrow = n, ncol = k)

  for (j in seq_len(k)) {
    diffs <- sweep(score_matrix, 2, centers_matrix[j, ], FUN = "-")
    dist_mat[, j] <- rowSums(diffs^2)
  }

  max.col(-dist_mat, ties.method = "first")
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
      warning("Could not read ", file_path, ": ", conditionMessage(e))
      return(NULL)
    }
  )

  empty_result <- list(
    metadata = tibble(),
    shape_matrix = matrix(numeric(0), nrow = 0, ncol = SHAPE_N_POINTS),
    display_matrix = matrix(numeric(0), nrow = 0, ncol = SHAPE_N_POINTS)
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

  missing_essential <- setdiff(essential, names(dat))
  if (length(missing_essential) > 0) {
    warning(
      "Skipping ", source_file,
      ". Missing columns: ",
      paste(missing_essential, collapse = ", ")
    )
    return(empty_result)
  }

  # Need either a directly stored I/I0 trace or current + baseline.
  has_iio_trace <- "I_over_I0_trace_pct" %in% names(dat)
  has_current <- all(c("current_pA", "baseline_pA") %in% names(dat))

  if (!has_iio_trace && !has_current) {
    warning(
      "Skipping ", source_file,
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
    sh$I_over_I0_trace_pct <- safe_numeric(dat$I_over_I0_trace_pct)
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

  # Convert possible 0-1 I/I0 values to percentage.
  q95_event <- suppressWarnings(stats::quantile(
    sh$event_I_over_I0_pct,
    0.95,
    na.rm = TRUE
  ))

  if (is.finite(q95_event) && q95_event <= 1.2) {
    sh$event_I_over_I0_pct <- sh$event_I_over_I0_pct * 100
  }

  q95_trace <- suppressWarnings(stats::quantile(
    sh$I_over_I0_trace_pct,
    0.95,
    na.rm = TRUE
  ))

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
    dplyr::filter(event_id %in% event_meta$event_id)

  # Use only the in-event part when available.
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

  event_splits <- split(sh_inside, sh_inside$event_id)
  grid <- seq(0, 1, length.out = SHAPE_N_POINTS)

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

    # Need at least 4 unique time points.
    if (length(unique(tr$time_relative_ms)) < 4) {
      next
    }

    x0 <- min(tr$time_relative_ms, na.rm = TRUE)
    x1 <- max(tr$time_relative_ms, na.rm = TRUE)

    if (!is.finite(x0) || !is.finite(x1) || x1 <= x0) {
      next
    }

    x_norm <- (tr$time_relative_ms - x0) / (x1 - x0)

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

    if (is.null(iio_interp) || any(!is.finite(iio_interp))) {
      next
    }

    # Shape signal for clustering:
    # 100-I/I0 converts the blockade into a positive deflection.
    # Divide by the event's own blockade amplitude so clustering is driven
    # mainly by waveform FORM rather than blockade depth.
    blockade <- pmax(0, 100 - iio_interp)

    amp <- suppressWarnings(stats::quantile(
      blockade,
      probs = 0.95,
      na.rm = TRUE,
      names = FALSE
    ))

    if (!is.finite(amp) || amp <= 1e-8) {
      next
    }

    shape_norm <- blockade / amp

    event_id_num <- safe_numeric(event_name)
    meta_row <- event_meta |>
      dplyr::filter(event_id == event_id_num) |>
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

  metadata <- dplyr::bind_rows(metadata_list)
  shape_matrix <- do.call(rbind, shape_list)
  display_matrix <- do.call(rbind, display_list)

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
  purrr::map(processed, "metadata")
)

shape_matrix_all <- do.call(
  rbind,
  purrr::map(processed, "shape_matrix")
)

display_matrix_all <- do.call(
  rbind,
  purrr::map(processed, "display_matrix")
)

if (nrow(metadata_all) < 10) {
  stop(
    "Too few valid filtered events for global clustering: ",
    nrow(metadata_all)
  )
}

if (
  nrow(shape_matrix_all) != nrow(metadata_all) ||
  nrow(display_matrix_all) != nrow(metadata_all)
) {
  stop("Internal row mismatch between metadata and waveform matrices.")
}

message("Valid filtered events: ", nrow(metadata_all))
message("Peptides: ", length(unique(metadata_all$peptide)))
message("Pores: ", nrow(metadata_all |> distinct(peptide, pore_id)))


# ============================================================
# 7. BUILD A BALANCED TRAINING SET
# ============================================================

training_meta <- metadata_all |>
  dplyr::group_by(peptide, pore_id) |>
  dplyr::group_modify(
    ~ dplyr::slice_sample(
      .x,
      n = min(nrow(.x), MAX_TRAIN_PER_PORE)
    )
  ) |>
  dplyr::ungroup()

training_meta <- training_meta |>
  dplyr::group_by(peptide) |>
  dplyr::group_modify(
    ~ dplyr::slice_sample(
      .x,
      n = min(nrow(.x), MAX_TRAIN_PER_PEPTIDE)
    )
  ) |>
  dplyr::ungroup()

if (nrow(training_meta) > MAX_TRAIN_TOTAL) {
  training_meta <- training_meta |>
    dplyr::slice_sample(n = MAX_TRAIN_TOTAL)
}

training_idx <- match(
  training_meta$event_uid,
  metadata_all$event_uid
)

training_shape_matrix <- shape_matrix_all[training_idx, , drop = FALSE]

message("Balanced training events: ", nrow(training_shape_matrix))


# ============================================================
# 8. PCA: REDUCE WAVEFORM DIMENSION BEFORE CLUSTERING
# ============================================================

pca_fit <- stats::prcomp(
  training_shape_matrix,
  center = TRUE,
  scale. = FALSE
)

var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
cum_var <- cumsum(var_explained)

n_pc <- which(cum_var >= PCA_VARIANCE_TARGET)[1]
if (is.na(n_pc)) {
  n_pc <- length(var_explained)
}

n_pc <- min(n_pc, PCA_MAX_COMPONENTS, ncol(pca_fit$x))
n_pc <- max(1, n_pc)

message(
  "PCA components retained: ", n_pc,
  " (cumulative variance = ",
  sprintf("%.1f%%", 100 * cum_var[n_pc]),
  ")"
)

training_scores <- pca_fit$x[, seq_len(n_pc), drop = FALSE]

pca_table <- tibble(
  PC = seq_along(var_explained),
  variance_explained = var_explained,
  cumulative_variance = cum_var
)

write.csv(
  pca_table,
  file.path(OUTPUT_FOLDER, "PCA_variance_explained.csv"),
  row.names = FALSE
)


# ============================================================
# 9. CHOOSE THE NUMBER OF GLOBAL CLUSTERS USING SILHOUETTE
# ============================================================

max_valid_k <- min(K_MAX, nrow(training_scores) - 1)
k_candidates <- seq.int(K_MIN, max_valid_k)

if (length(k_candidates) == 0) {
  stop("Not enough training events to test clustering.")
}

training_dist <- stats::dist(training_scores)

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
      return(tibble(k = k, average_silhouette = NA_real_))
    }

    sil <- cluster::silhouette(
      km_try$cluster,
      training_dist
    )

    tibble(
      k = k,
      average_silhouette = mean(sil[, "sil_width"], na.rm = TRUE)
    )
  }
)

valid_k_results <- k_results |>
  dplyr::filter(is.finite(average_silhouette))

if (nrow(valid_k_results) == 0) {
  stop("Could not obtain a valid silhouette score for any tested k.")
}

best_k <- valid_k_results |>
  dplyr::arrange(dplyr::desc(average_silhouette), k) |>
  dplyr::slice(1) |>
  dplyr::pull(k)

message("Selected global number of clusters: k = ", best_k)
print(k_results)

write.csv(
  k_results,
  file.path(OUTPUT_FOLDER, "cluster_number_silhouette_scores.csv"),
  row.names = FALSE
)

p_k <- ggplot2::ggplot(
  k_results,
  ggplot2::aes(x = k, y = average_silhouette)
) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_vline(
    xintercept = best_k,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::scale_x_continuous(breaks = k_candidates) +
  ggplot2::labs(
    x = "Number of clusters (k)",
    y = "Average silhouette score",
    title = paste0("Global event-shape clustering: selected k = ", best_k)
  ) +
  ggplot2::theme_classic(base_size = 10)

ggplot2::ggsave(
  file.path(OUTPUT_FOLDER, "01_cluster_number_selection.png"),
  p_k,
  width = 6.0,
  height = 4.2,
  dpi = PNG_DPI,
  bg = "white"
)


# ============================================================
# 10. FIT FINAL GLOBAL K-MEANS MODEL
# ============================================================

final_kmeans <- stats::kmeans(
  training_scores,
  centers = best_k,
  nstart = 100,
  iter.max = 300
)

# Project ALL filtered waveforms into the training PCA space.
all_scores_full <- stats::predict(
  pca_fit,
  newdata = shape_matrix_all
)

all_scores <- all_scores_full[, seq_len(n_pc), drop = FALSE]

raw_cluster <- nearest_centroid_assignment(
  score_matrix = all_scores,
  centers_matrix = final_kmeans$centers
)

metadata_all$cluster_raw <- raw_cluster

# Relabel clusters by global abundance so Cluster 1 is the most common
# cluster across all filtered events.
cluster_order <- metadata_all |>
  dplyr::count(cluster_raw, name = "N") |>
  dplyr::arrange(dplyr::desc(N), cluster_raw) |>
  dplyr::mutate(cluster = dplyr::row_number())

relabel_map <- setNames(
  cluster_order$cluster,
  cluster_order$cluster_raw
)

metadata_all$cluster <- as.integer(
  relabel_map[as.character(metadata_all$cluster_raw)]
)

metadata_all$cluster_label <- factor(
  paste0("Cluster ", metadata_all$cluster),
  levels = paste0("Cluster ", seq_len(best_k))
)


# ============================================================
# 11. EVENT-LEVEL ASSIGNMENT TABLE
# ============================================================

assignment_file <- file.path(
  OUTPUT_FOLDER,
  "all_filtered_events_global_cluster_assignments.csv"
)

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
  assignment_file,
  row.names = FALSE
)


# ============================================================
# 12. CLUSTER SUMMARY + REPRESENTATIVE REAL EVENT
# ============================================================

# Distance of each event to the centroid of its assigned RAW cluster.
distance_to_centroid <- rep(NA_real_, nrow(all_scores))

for (j in seq_len(nrow(final_kmeans$centers))) {
  idx <- which(metadata_all$cluster_raw == j)
  if (length(idx) == 0) next

  diffs <- sweep(
    all_scores[idx, , drop = FALSE],
    2,
    final_kmeans$centers[j, ],
    FUN = "-"
  )

  distance_to_centroid[idx] <- sqrt(rowSums(diffs^2))
}

metadata_all$distance_to_cluster_centroid <- distance_to_centroid

representative_events <- metadata_all |>
  dplyr::group_by(cluster, cluster_label) |>
  dplyr::arrange(distance_to_cluster_centroid, .by_group = TRUE) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

write.csv(
  representative_events,
  file.path(OUTPUT_FOLDER, "representative_real_event_per_cluster.csv"),
  row.names = FALSE
)

cluster_summary <- metadata_all |>
  dplyr::group_by(cluster, cluster_label) |>
  dplyr::summarise(
    N_events = dplyr::n(),
    N_peptides = dplyr::n_distinct(peptide),
    N_pores = dplyr::n_distinct(paste(peptide, pore_id)),
    median_dwell_ms = stats::median(dwell_total_ms, na.rm = TRUE),
    median_Ires_pct = stats::median(Ires_pct, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  cluster_summary,
  file.path(OUTPUT_FOLDER, "global_cluster_summary.csv"),
  row.names = FALSE
)


# ============================================================
# 13. REPRESENTATIVE SHAPE OF EACH GLOBAL CLUSTER
# ============================================================

norm_time <- seq(0, 1, length.out = SHAPE_N_POINTS)

# Long version of actual I/I0 traces used only for visualisation.
display_long <- tibble(
  row_index = rep(seq_len(nrow(display_matrix_all)), each = SHAPE_N_POINTS),
  normalised_time = rep(norm_time, times = nrow(display_matrix_all)),
  I_over_I0_pct = as.vector(t(display_matrix_all))
) |>
  dplyr::mutate(
    event_uid = metadata_all$event_uid[row_index],
    cluster = metadata_all$cluster[row_index],
    cluster_label = metadata_all$cluster_label[row_index]
  )

cluster_median_shapes <- display_long |>
  dplyr::group_by(cluster, cluster_label, normalised_time) |>
  dplyr::summarise(
    median_I_over_I0_pct = stats::median(I_over_I0_pct, na.rm = TRUE),
    Q25 = stats::quantile(I_over_I0_pct, 0.25, na.rm = TRUE),
    Q75 = stats::quantile(I_over_I0_pct, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  cluster_median_shapes,
  file.path(OUTPUT_FOLDER, "global_cluster_median_shapes.csv"),
  row.names = FALSE
)

# Select real traces nearest to the cluster centroid for grey overlays.
overlay_event_ids <- metadata_all |>
  dplyr::group_by(cluster, cluster_label) |>
  dplyr::arrange(distance_to_cluster_centroid, .by_group = TRUE) |>
  dplyr::slice_head(n = MAX_OVERLAY_PER_CLUSTER) |>
  dplyr::ungroup() |>
  dplyr::pull(event_uid)

overlay_long <- display_long |>
  dplyr::filter(event_uid %in% overlay_event_ids)

p_shapes <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = overlay_long,
    ggplot2::aes(
      x = normalised_time,
      y = I_over_I0_pct,
      group = event_uid
    ),
    colour = "grey78",
    linewidth = 0.30,
    alpha = 0.65
  ) +
  ggplot2::geom_line(
    data = cluster_median_shapes,
    ggplot2::aes(
      x = normalised_time,
      y = median_I_over_I0_pct
    ),
    colour = "black",
    linewidth = 0.85
  ) +
  ggplot2::facet_wrap(~ cluster_label, nrow = 1, scales = "free_y") +
  ggplot2::labs(
    x = "Normalised event time",
    y = expression(I/I[0] ~ " (%)"),
    title = "Representative waveform of each global event-shape cluster"
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "none"
  )

ggplot2::ggsave(
  file.path(OUTPUT_FOLDER, "02_global_cluster_representative_shapes.png"),
  p_shapes,
  width = max(8, 3.0 * best_k),
  height = 4.0,
  dpi = PNG_DPI,
  bg = "white",
  limitsize = FALSE
)


# ============================================================
# 14. IRES VS DWELL SCATTER COLOURED BY GLOBAL CLUSTER
# ============================================================

p_scatter <- ggplot2::ggplot(
  metadata_all,
  ggplot2::aes(
    x = Ires_pct,
    y = dwell_total_ms,
    colour = cluster_label
  )
) +
  ggplot2::geom_point(
    alpha = 0.50,
    size = 0.65
  ) +
  ggplot2::scale_y_log10() +
  ggplot2::facet_wrap(~ peptide, scales = "free") +
  ggplot2::labs(
    x = expression(I/I[0] ~ " (%)"),
    y = "Dwell time (ms)",
    colour = "Shape cluster",
    title = "Global waveform clusters mapped onto I/I0 vs dwell time"
  ) +
  ggplot2::theme_classic(base_size = 8) +
  ggplot2::theme(
    legend.position = "bottom",
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(OUTPUT_FOLDER, "03_Ires_vs_dwell_by_global_cluster.png"),
  p_scatter,
  width = 12,
  height = max(5.5, 2.6 * ceiling(length(unique(metadata_all$peptide)) / 3)),
  dpi = PNG_DPI,
  bg = "white",
  limitsize = FALSE
)


# ============================================================
# 15. CLUSTER PROPORTIONS PER PORE
# ============================================================

all_cluster_levels <- paste0("Cluster ", seq_len(best_k))

pore_props <- metadata_all |>
  dplyr::count(peptide, pore_id, cluster_label, name = "N_events") |>
  tidyr::complete(
    peptide,
    pore_id,
    cluster_label = factor(
      all_cluster_levels,
      levels = all_cluster_levels
    ),
    fill = list(N_events = 0)
  ) |>
  dplyr::group_by(peptide, pore_id) |>
  dplyr::mutate(
    N_total_filtered = sum(N_events),
    proportion = dplyr::if_else(
      N_total_filtered > 0,
      N_events / N_total_filtered,
      0
    )
  ) |>
  dplyr::ungroup()

write.csv(
  pore_props,
  file.path(OUTPUT_FOLDER, "cluster_proportions_by_pore.csv"),
  row.names = FALSE
)


# ============================================================
# 16. PEPTIDE COMPARISON - TWO VERSIONS
# ============================================================

# A) Event-weighted: every event has equal weight.
peptide_props_event_weighted <- metadata_all |>
  dplyr::count(peptide, cluster_label, name = "N_events") |>
  tidyr::complete(
    peptide,
    cluster_label = factor(
      all_cluster_levels,
      levels = all_cluster_levels
    ),
    fill = list(N_events = 0)
  ) |>
  dplyr::group_by(peptide) |>
  dplyr::mutate(
    N_total_filtered = sum(N_events),
    proportion = N_events / N_total_filtered
  ) |>
  dplyr::ungroup()

write.csv(
  peptide_props_event_weighted,
  file.path(OUTPUT_FOLDER, "cluster_proportions_by_peptide_event_weighted.csv"),
  row.names = FALSE
)

# B) Pore-balanced: calculate cluster percentages in each pore first,
# then average those percentages across pores of the same peptide.
# This prevents one pore with many more events from dominating the peptide.
peptide_props_pore_balanced <- pore_props |>
  dplyr::group_by(peptide, cluster_label) |>
  dplyr::summarise(
    mean_proportion = mean(proportion, na.rm = TRUE),
    sd_proportion = stats::sd(proportion, na.rm = TRUE),
    N_pores = dplyr::n_distinct(pore_id),
    .groups = "drop"
  )

write.csv(
  peptide_props_pore_balanced,
  file.path(OUTPUT_FOLDER, "cluster_proportions_by_peptide_pore_balanced.csv"),
  row.names = FALSE
)


# ============================================================
# 17. FINAL PEPTIDE COMPARISON FIGURE
# ============================================================

p_prop <- ggplot2::ggplot(
  peptide_props_pore_balanced,
  ggplot2::aes(
    x = peptide,
    y = mean_proportion,
    fill = cluster_label
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.20
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Mean proportion of filtered events",
    fill = "Shape cluster",
    title = "Global event-shape cluster distribution across peptides",
    subtitle = "Pore-balanced proportions"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(OUTPUT_FOLDER, "04_global_cluster_proportions_across_peptides.png"),
  p_prop,
  width = max(8, 0.80 * length(unique(metadata_all$peptide)) + 4),
  height = 5.2,
  dpi = PNG_DPI,
  bg = "white",
  limitsize = FALSE
)


# ============================================================
# 18. SMALL TEXT SUMMARY
# ============================================================

summary_lines <- c(
  paste0("Root folder: ", ROOT_FOLDER),
  paste0("Number of event-shape CSV files: ", length(shape_files)),
  paste0("Number of peptides: ", length(unique(metadata_all$peptide))),
  paste0("Number of peptide/pore combinations: ", nrow(metadata_all |> distinct(peptide, pore_id))),
  paste0("Filtered events used for assignment: ", nrow(metadata_all)),
  paste0("Filter: dwell >= ", MIN_DWELL_MS, " ms and Ires = ", MIN_IRES_PCT, "-", MAX_IRES_PCT, "%"),
  paste0("Balanced events used to train clustering: ", nrow(training_meta)),
  paste0("PCA components retained: ", n_pc),
  paste0("Selected number of global clusters: ", best_k),
  "",
  "Interpretation:",
  "- Global cluster labels have the same meaning across all peptides.",
  "- Clustering uses normalised waveform shape, not raw blockade amplitude.",
  "- Representative cluster traces are extracted after clustering.",
  "- The final peptide comparison uses pore-balanced cluster proportions."
)

writeLines(
  summary_lines,
  con = file.path(OUTPUT_FOLDER, "README_global_shape_clustering.txt")
)


# ============================================================
# 19. FINISHED
# ============================================================

message("")
message("============================================================")
message("GLOBAL EVENT-SHAPE CLUSTERING COMPLETE")
message("============================================================")
message("Selected k = ", best_k)
message("All outputs saved in:")
message(OUTPUT_FOLDER)
message("")
message("Main figures:")
message("  01_cluster_number_selection.png")
message("  02_global_cluster_representative_shapes.png")
message("  03_Ires_vs_dwell_by_global_cluster.png")
message("  04_global_cluster_proportions_across_peptides.png")
message("")
message("Main tables:")
message("  all_filtered_events_global_cluster_assignments.csv")
message("  representative_real_event_per_cluster.csv")
message("  global_cluster_summary.csv")
message("  cluster_proportions_by_pore.csv")
message("  cluster_proportions_by_peptide_event_weighted.csv")
message("  cluster_proportions_by_peptide_pore_balanced.csv")
message("============================================================")
