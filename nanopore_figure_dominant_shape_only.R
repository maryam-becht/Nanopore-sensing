# ============================================================
# ONE PEPTIDE FOLDER AT A TIME
# Final version: panel e = dominant event shape
#
# What panel e shows:
#   - filter events with dwell >= 0.5 ms and Ires <= 85%
#   - compare waveform SHAPES across filtered events
#   - identify the most frequent shape family
#   - display a light-grey overlay of that dominant family
#   - display the median dominant shape in black
#
# Expected structure in one peptide folder:
#   peptide_name/
#     peptide_P1_nanopore_R_data.xlsx
#     peptide_P1_nanopore_R_data_event_shapes.csv
#     peptide_P2_nanopore_R_data.xlsx
#     peptide_P2_nanopore_R_data_event_shapes.csv
#     ...
#
# Output:
#   peptide_name_R_folder_output/
#     peptide_name_dominant_shape_only.png
#     peptide_name_dominant_shape_summary.csv
# ============================================================

# ============================================================
# 0. PACKAGES
# ============================================================

required_packages <- c(
  "openxlsx",
  "dplyr",
  "purrr",
  "ggplot2",
  "patchwork",
  "tibble",
  "grid"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

library(openxlsx)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(tibble)
library(grid)


# ============================================================
# 1. SELECT ONE PEPTIDE FOLDER
# ============================================================

if (
  interactive() &&
  .Platform$OS.type == "windows"
) {

  peptide_folder <- choose.dir(
    default = "C:/Users/mbech/Documents/Peptides",
    caption = "Choisis le dossier d'un peptide"
  )

} else {

  peptide_folder <- readline(
    "Entre le chemin complet du dossier peptide : "
  )
}

if (
  is.na(peptide_folder) ||
  peptide_folder == ""
) {
  stop("Aucun dossier sélectionné.")
}

peptide_folder <- normalizePath(
  peptide_folder,
  winslash = "/",
  mustWork = TRUE
)

peptide_id <- basename(peptide_folder)

output_folder <- file.path(
  peptide_folder,
  paste0(peptide_id, "_R_folder_output")
)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

message("Dossier peptide : ", peptide_folder)
message("Peptide : ", peptide_id)


# ============================================================
# 2. PARAMETERS
# ============================================================

MIN_DWELL_MS <- 0.2
MAX_DWELL_MS <- 100
MIN_I_OVER_I0 <- 0
MAX_I_OVER_I0 <- 100
DISPLAY_MAX_DWELL_HIST <- 5
MIN_EVENTS_PER_PORE <- 100

DWELL_BINWIDTH <- 0.20
IIO_BINWIDTH <- 2

SCATTER_MAX_POINTS <- 5000
TRACE_MAX_POINTS <- 4500

# Dominant-shape filters
DOM_MIN_DWELL_MS <- 0.5
DOM_MAX_IRES_PCT <- 85

# Shape-comparison parameters
SHAPE_COMPARE_NPTS <- 80
SHAPE_K_MAX <- 3
DOM_OVERLAY_MAX_TRACES <- 18
DOM_DISPLAY_NPTS <- 180

# Dominant-shape display
DOM_TIME_SCALEBAR_MS <- 0.5
DOM_CURRENT_SCALEBAR_PA <- 10
DOM_LEFT_PAD_MS <- 0.60
DOM_RIGHT_PAD_MS <- 0.90
DOM_TITLE_SIZE <- 3.0

PNG_DPI <- 180
FIGURE_WIDTH_IN <- 16.0
HEIGHT_PER_PORE_IN <- 2.15

set.seed(123)


# ============================================================
# 3. FIND ALL EXCEL FILES = ALL PORES OF THIS PEPTIDE
# ============================================================

excel_files <- list.files(
  path = peptide_folder,
  pattern = "\\.(xlsx|xlsm)$",
  full.names = TRUE,
  ignore.case = TRUE
)

excel_files <- excel_files[
  !grepl("^~\\$", basename(excel_files))
]

excel_files <- sort(excel_files)

if (length(excel_files) == 0) {
  stop(
    "Aucun fichier Excel trouvé dans le dossier : ",
    peptide_folder
  )
}

message("Nombre de fichiers pore trouvés : ", length(excel_files))
print(basename(excel_files))


# ============================================================
# 4. HELPERS
# ============================================================

extract_pore_id_from_filename <- function(file_path) {

  fname <- tools::file_path_sans_ext(basename(file_path))

  pore_id <- sub(
    ".*_(P[0-9]+)_nanopore_R_data$",
    "\\1",
    fname,
    ignore.case = TRUE
  )

  if (
    identical(pore_id, fname) ||
    is.na(pore_id) ||
    pore_id == ""
  ) {

    m <- regexpr("P[0-9]+", fname, ignore.case = TRUE)

    if (m[1] > 0) {
      pore_id <- regmatches(fname, m)
    } else {
      pore_id <- fname
    }
  }

  toupper(pore_id)
}

clean_names_simple <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

match_first_column <- function(df, candidates) {
  nm <- clean_names_simple(names(df))
  idx <- match(clean_names_simple(candidates), nm)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    return(NA_integer_)
  }
  idx[1]
}

find_sheet_name <- function(wb, candidates, pattern = NULL) {

  sheet_names <- names(wb)
  sheet_names_clean <- clean_names_simple(sheet_names)

  candidate_clean <- clean_names_simple(candidates)
  hit <- match(candidate_clean, sheet_names_clean)
  hit <- hit[!is.na(hit)]

  if (length(hit) >= 1) {
    return(sheet_names[hit[1]])
  }

  if (!is.null(pattern)) {
    hit2 <- grep(pattern, sheet_names, ignore.case = TRUE)
    if (length(hit2) >= 1) {
      return(sheet_names[hit2[1]])
    }
  }

  NA_character_
}

scale_density_to_counts <- function(x, binwidth, x_from, x_to, n_events) {

  if (length(x) < 10) {
    return(tibble(x = numeric(0), scaled_count = numeric(0)))
  }

  d <- density(
    x,
    from = x_from,
    to = x_to,
    n = 300,
    na.rm = TRUE
  )

  tibble(
    x = d$x,
    scaled_count = d$y * n_events * binwidth
  )
}

fit_shifted_exponential <- function(x_ms, threshold_ms) {

  x_ms <- x_ms[is.finite(x_ms)]
  x_ms <- x_ms[x_ms >= threshold_ms]

  shifted <- x_ms - threshold_ms

  if (length(shifted) < 5 || mean(shifted) <= 0) {
    return(tibble(x = numeric(0), scaled_count = numeric(0)))
  }

  lambda_hat <- 1 / mean(shifted)

  x_grid <- seq(
    threshold_ms,
    DISPLAY_MAX_DWELL_HIST,
    length.out = 250
  )

  density_fit <- lambda_hat * exp(-lambda_hat * (x_grid - threshold_ms))

  tibble(
    x = x_grid,
    scaled_count = density_fit * length(x_ms) * DWELL_BINWIDTH
  )
}

moving_average3 <- function(x) {
  stats::filter(x, rep(1/3, 3), sides = 2, circular = FALSE) |>
    as.numeric()
}


# ============================================================
# 5. READ EVENTS AND TRACE DATA
# ============================================================

read_events_one_file <- function(file_path, pore_name) {

  wb <- openxlsx::loadWorkbook(file_path)

  events_sheet <- find_sheet_name(
    wb,
    candidates = c("events", "events_filtered"),
    pattern = "events"
  )

  if (is.na(events_sheet)) {
    stop("Aucune feuille events trouvée dans ", basename(file_path))
  }

  ev <- openxlsx::read.xlsx(
    wb,
    sheet = events_sheet,
    check.names = FALSE
  )

  if (nrow(ev) == 0) {
    return(tibble(
      pore_id = character(0),
      dwell_total_ms = numeric(0),
      I_over_I0_pct = numeric(0)
    ))
  }

  dwell_idx <- match_first_column(
    ev,
    c(
      "dwell_total_ms",
      "dwell_ms",
      "duration_ms",
      "event_duration_ms",
      "dwell_time_ms",
      "duration_s",
      "duration"
    )
  )

  iio_idx <- match_first_column(
    ev,
    c(
      "I_over_I0_pct",
      "i_over_i0_pct",
      "i_i0_pct",
      "event_I_over_I0_pct",
      "ires_pct",
      "residual_current_pct",
      "I_over_I0",
      "i_over_i0"
    )
  )

  if (is.na(dwell_idx) || is.na(iio_idx)) {
    stop(
      "Colonnes dwell/I_over_I0 manquantes dans ",
      basename(file_path),
      " / ",
      events_sheet
    )
  }

  dwell <- as.numeric(ev[[dwell_idx]])
  iio <- as.numeric(ev[[iio_idx]])

  # If dwell column is in seconds, convert to ms.
  if (grepl("_s$|^duration$", clean_names_simple(names(ev)[dwell_idx]))) {
    dwell <- dwell * 1000
  } else if (median(dwell, na.rm = TRUE) < 0.05) {
    # safety net: values probably in seconds
    dwell <- dwell * 1000
  }

  # If I/I0 is 0-1, convert to percent.
  if (stats::quantile(iio, 0.95, na.rm = TRUE) <= 1.2) {
    iio <- iio * 100
  }

  tibble(
    pore_id = pore_name,
    dwell_total_ms = dwell,
    I_over_I0_pct = iio
  ) |>
    dplyr::filter(
      is.finite(dwell_total_ms),
      is.finite(I_over_I0_pct)
    )
}

read_trace_one_file <- function(file_path, pore_name) {

  wb <- openxlsx::loadWorkbook(file_path)

  trace_sheet <- find_sheet_name(
    wb,
    candidates = c(paste0("trace_", pore_name), "trace"),
    pattern = "trace"
  )

  if (is.na(trace_sheet)) {
    stop("Aucune feuille trace trouvée dans ", basename(file_path))
  }

  tr <- openxlsx::read.xlsx(
    wb,
    sheet = trace_sheet,
    check.names = FALSE
  )

  if (nrow(tr) == 0) {
    return(tibble(
      pore_id = character(0),
      time_s = numeric(0),
      current_pA = numeric(0),
      baseline_pA = numeric(0)
    ))
  }

  time_idx <- match_first_column(tr, c("time_s", "time_sec", "time"))
  current_idx <- match_first_column(tr, c("current_pA", "current_pa", "current"))
  baseline_idx <- match_first_column(tr, c("baseline_pA", "baseline_pa", "baseline"))

  if (is.na(time_idx) || is.na(current_idx)) {
    stop(
      "Colonnes time/current manquantes dans ",
      basename(file_path),
      " / ",
      trace_sheet
    )
  }

  out <- tibble(
    pore_id = pore_name,
    time_s = as.numeric(tr[[time_idx]]),
    current_pA = as.numeric(tr[[current_idx]])
  ) |>
    dplyr::filter(is.finite(time_s), is.finite(current_pA))

  baseline_value <- if (!is.na(baseline_idx)) {
    stats::median(as.numeric(tr[[baseline_idx]]), na.rm = TRUE)
  } else {
    stats::median(out$current_pA, na.rm = TRUE)
  }

  out <- out |>
    dplyr::mutate(baseline_pA = baseline_value)

  if (nrow(out) > TRACE_MAX_POINTS) {
    keep_idx <- unique(round(seq(1, nrow(out), length.out = TRACE_MAX_POINTS)))
    out <- out[keep_idx, , drop = FALSE]
  }

  out
}


# ============================================================
# 6. DOMINANT SHAPE READER
# ============================================================

read_dominant_shape_one_file <- function(file_path, pore_name) {

  shapes_file <- paste0(
    tools::file_path_sans_ext(file_path),
    "_event_shapes.csv"
  )

  empty_result <- list(
    overlay = tibble(
      event_id = integer(0),
      time_relative_ms = numeric(0),
      current_centered_pA = numeric(0)
    ),
    median_trace = tibble(
      time_relative_ms = numeric(0),
      current_centered_pA = numeric(0)
    ),
    summary = tibble(
      pore_id = pore_name,
      n_filtered = NA_integer_,
      n_in_dominant_cluster = NA_integer_,
      k_used = NA_integer_,
      status = NA_character_
    )
  )

  if (!file.exists(shapes_file)) {
    empty_result$summary$status <- "event_shapes CSV missing"
    return(empty_result)
  }

  shapes <- tryCatch(
    utils::read.csv(shapes_file, check.names = FALSE),
    error = function(e) NULL
  )

  if (is.null(shapes)) {
    empty_result$summary$status <- "could not read event_shapes CSV"
    return(empty_result)
  }

  required_shape_columns <- c(
    "event_id",
    "time_relative_ms",
    "current_pA",
    "baseline_pA",
    "dwell_total_ms",
    "event_I_over_I0_pct"
  )

  missing_cols <- setdiff(required_shape_columns, names(shapes))
  if (length(missing_cols) > 0) {
    empty_result$summary$status <- paste0(
      "missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
    return(empty_result)
  }

  sh <- shapes |>
    dplyr::transmute(
      event_id = as.integer(event_id),
      time_relative_ms = as.numeric(time_relative_ms),
      current_pA = as.numeric(current_pA),
      baseline_pA = as.numeric(baseline_pA),
      dwell_total_ms = as.numeric(dwell_total_ms),
      event_I_over_I0_pct = as.numeric(event_I_over_I0_pct)
    ) |>
    dplyr::filter(
      is.finite(event_id),
      is.finite(time_relative_ms),
      is.finite(current_pA),
      is.finite(baseline_pA),
      is.finite(dwell_total_ms),
      is.finite(event_I_over_I0_pct)
    ) |>
    dplyr::mutate(
      current_centered_pA = current_pA - baseline_pA
    )

  event_summary <- sh |>
    dplyr::distinct(event_id, dwell_total_ms, event_I_over_I0_pct) |>
    dplyr::filter(
      dwell_total_ms >= DOM_MIN_DWELL_MS,
      dwell_total_ms <= MAX_DWELL_MS,
      event_I_over_I0_pct >= MIN_I_OVER_I0,
      event_I_over_I0_pct <= DOM_MAX_IRES_PCT
    )

  n_filtered <- nrow(event_summary)

  if (n_filtered < 2) {
    empty_result$summary$n_filtered <- n_filtered
    empty_result$summary$status <- "not enough filtered events"
    return(empty_result)
  }

  # Build shape profiles for comparison only.
  compare_grid <- seq(0, 1, length.out = SHAPE_COMPARE_NPTS)

  shape_profiles <- purrr::map(
    event_summary$event_id,
    function(this_event_id) {

      tr <- sh |>
        dplyr::filter(event_id == this_event_id) |>
        dplyr::arrange(time_relative_ms)

      if (nrow(tr) < 6) {
        return(NULL)
      }

      x <- tr$time_relative_ms
      y <- tr$current_centered_pA

      x0 <- min(x, na.rm = TRUE)
      x1 <- max(x, na.rm = TRUE)

      if (!is.finite(x0) || !is.finite(x1) || x1 <= x0) {
        return(NULL)
      }

      x_scaled <- (x - x0) / (x1 - x0)

      # Normalize amplitude for SHAPE comparison only.
      amp <- abs(min(y, na.rm = TRUE))
      if (!is.finite(amp) || amp <= 0) {
        return(NULL)
      }
      y_norm <- y / amp

      interp <- approx(
        x = x_scaled,
        y = y_norm,
        xout = compare_grid,
        ties = mean,
        rule = 2
      )$y

      interp
    }
  )

  valid_idx <- which(vapply(shape_profiles, Negate(is.null), logical(1)))

  if (length(valid_idx) < 2) {
    empty_result$summary$n_filtered <- n_filtered
    empty_result$summary$status <- "not enough valid shape profiles"
    return(empty_result)
  }

  valid_event_ids <- event_summary$event_id[valid_idx]
  shape_matrix <- do.call(rbind, shape_profiles[valid_idx])
  rownames(shape_matrix) <- valid_event_ids

  if (nrow(shape_matrix) == 2) {
    cluster_id <- c(1, 1)
    k_used <- 1
  } else {
    k_used <- min(SHAPE_K_MAX, nrow(shape_matrix))
    km <- stats::kmeans(shape_matrix, centers = k_used, nstart = 20)
    cluster_id <- km$cluster
  }

  cluster_tbl <- tibble(
    event_id = as.integer(rownames(shape_matrix)),
    cluster = cluster_id
  ) |>
    dplyr::count(cluster, name = "cluster_n") |>
    dplyr::arrange(dplyr::desc(cluster_n), cluster) |>
    dplyr::slice(1)

  dominant_cluster <- cluster_tbl$cluster[1]
  dominant_n <- cluster_tbl$cluster_n[1]

  dominant_event_ids <- tibble(
    event_id = as.integer(rownames(shape_matrix)),
    cluster = cluster_id
  ) |>
    dplyr::filter(cluster == dominant_cluster) |>
    dplyr::pull(event_id)

  dominant_shape_matrix <- shape_matrix[as.character(dominant_event_ids), , drop = FALSE]

  median_profile <- apply(dominant_shape_matrix, 2, stats::median, na.rm = TRUE)

  profile_distances <- apply(
    dominant_shape_matrix,
    1,
    function(z) sqrt(mean((z - median_profile)^2, na.rm = TRUE))
  )

  selected_overlay_ids <- names(sort(profile_distances))[seq_len(min(DOM_OVERLAY_MAX_TRACES, length(profile_distances)))] |>
    as.integer()

  overlay_raw <- sh |>
    dplyr::filter(event_id %in% selected_overlay_ids) |>
    dplyr::arrange(event_id, time_relative_ms)

  dominant_raw <- sh |>
    dplyr::filter(event_id %in% dominant_event_ids) |>
    dplyr::arrange(event_id, time_relative_ms)

  x_left <- min(dominant_raw$time_relative_ms, na.rm = TRUE)
  x_right <- max(dominant_raw$time_relative_ms, na.rm = TRUE)

  x_left <- min(x_left, -DOM_LEFT_PAD_MS)
  x_right <- max(x_right, max(dominant_raw$dwell_total_ms, na.rm = TRUE) + DOM_RIGHT_PAD_MS)

  display_grid <- seq(x_left, x_right, length.out = DOM_DISPLAY_NPTS)

  display_profiles <- purrr::map_dfr(
    unique(dominant_raw$event_id),
    function(this_event_id) {
      tr <- dominant_raw |>
        dplyr::filter(event_id == this_event_id) |>
        dplyr::arrange(time_relative_ms)

      interp <- approx(
        x = tr$time_relative_ms,
        y = tr$current_centered_pA,
        xout = display_grid,
        rule = 2,
        ties = mean
      )$y

      tibble(
        event_id = this_event_id,
        time_relative_ms = display_grid,
        current_centered_pA = interp
      )
    }
  )

  median_trace <- display_profiles |>
    dplyr::group_by(time_relative_ms) |>
    dplyr::summarise(
      current_centered_pA = stats::median(current_centered_pA, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      current_centered_pA = moving_average3(current_centered_pA)
    )

  overlay_trace <- overlay_raw |>
    dplyr::mutate(
      current_centered_pA = moving_average3(current_centered_pA)
    )

  summary_tbl <- tibble(
    pore_id = pore_name,
    n_filtered = n_filtered,
    n_in_dominant_cluster = dominant_n,
    k_used = k_used,
    status = "OK"
  )

  list(
    overlay = overlay_trace,
    median_trace = median_trace,
    summary = summary_tbl
  )
}


# ============================================================
# 7. READ ALL DATA
# ============================================================

pore_ids <- vapply(excel_files, extract_pore_id_from_filename, character(1))

pore_info <- tibble(
  file_path = excel_files,
  file_name = basename(excel_files),
  pore_id = pore_ids
) |>
  dplyr::mutate(
    pore_number = suppressWarnings(as.integer(sub("^P", "", pore_id)))
  ) |>
  dplyr::arrange(pore_number, pore_id) |>
  dplyr::distinct(pore_id, .keep_all = TRUE)

pore_names <- pore_info$pore_id
n_pores <- nrow(pore_info)

message("Nombre de pores détectés : ", n_pores)
print(pore_info |> dplyr::select(pore_id, file_name))

events_all <- purrr::map2_dfr(
  pore_info$file_path,
  pore_info$pore_id,
  read_events_one_file
)

trace_all <- purrr::map2_dfr(
  pore_info$file_path,
  pore_info$pore_id,
  read_trace_one_file
)

dominant_shape_results <- setNames(
  purrr::map2(
    pore_info$file_path,
    pore_info$pore_id,
    read_dominant_shape_one_file
  ),
  pore_info$pore_id
)

dominant_shape_summary <- purrr::map_dfr(dominant_shape_results, "summary")

write.csv(
  dominant_shape_summary,
  file.path(output_folder, paste0(peptide_id, "_dominant_shape_summary.csv")),
  row.names = FALSE
)


# ============================================================
# 8. FILTER EVENTS FOR FIGURE a-c
# ============================================================

events_plot <- events_all |>
  dplyr::filter(
    dwell_total_ms >= MIN_DWELL_MS,
    dwell_total_ms <= MAX_DWELL_MS,
    I_over_I0_pct >= MIN_I_OVER_I0,
    I_over_I0_pct <= MAX_I_OVER_I0
  )

pore_counts <- events_plot |>
  dplyr::count(pore_id, name = "n_events")

pores_keep <- pore_counts |>
  dplyr::filter(n_events >= MIN_EVENTS_PER_PORE) |>
  dplyr::pull(pore_id)

if (length(pores_keep) == 0) {
  stop("Aucun pore n'a au moins ", MIN_EVENTS_PER_PORE, " events.")
}

events_plot <- events_plot |>
  dplyr::filter(pore_id %in% pores_keep)

trace_all <- trace_all |>
  dplyr::filter(pore_id %in% pores_keep)

pore_info <- pore_info |>
  dplyr::filter(pore_id %in% pores_keep)

pore_names <- pore_info$pore_id
n_pores <- length(pore_names)

message("Pores conservés pour la figure : ", paste(pore_names, collapse = ", "))


# ============================================================
# 9. COMMON AXES FOR HISTOGRAMS
# ============================================================

dwell_hist_ymax <- max(
  purrr::map_dbl(
    pore_names,
    function(current_pore) {
      x <- events_plot |>
        dplyr::filter(
          pore_id == current_pore,
          dwell_total_ms >= MIN_DWELL_MS,
          dwell_total_ms <= DISPLAY_MAX_DWELL_HIST
        ) |>
        dplyr::pull(dwell_total_ms)

      if (length(x) == 0) return(0)

      breaks <- seq(MIN_DWELL_MS, DISPLAY_MAX_DWELL_HIST + DWELL_BINWIDTH, by = DWELL_BINWIDTH)
      h <- hist(x, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)
      curve <- fit_shifted_exponential(x, MIN_DWELL_MS)

      max(c(h$counts, curve$scaled_count), na.rm = TRUE)
    }
  ),
  na.rm = TRUE
) * 1.05

current_hist_ymax <- max(
  purrr::map_dbl(
    pore_names,
    function(current_pore) {
      x <- events_plot |>
        dplyr::filter(pore_id == current_pore) |>
        dplyr::pull(I_over_I0_pct)

      if (length(x) == 0) return(0)

      breaks <- seq(MIN_I_OVER_I0, MAX_I_OVER_I0 + IIO_BINWIDTH, by = IIO_BINWIDTH)
      h <- hist(x, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)
      curve <- scale_density_to_counts(x, IIO_BINWIDTH, MIN_I_OVER_I0, MAX_I_OVER_I0, length(x))

      max(c(h$counts, curve$scaled_count), na.rm = TRUE)
    }
  ),
  na.rm = TRUE
) * 1.05


# ============================================================
# 10. THEME
# ============================================================

article_theme <- ggplot2::theme_classic(base_size = 8) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 8, face = "plain", hjust = 0.5),
    axis.title = ggplot2::element_text(size = 8, colour = "black"),
    axis.text = ggplot2::element_text(size = 7, colour = "black"),
    plot.tag = ggplot2::element_text(size = 11, face = "bold", colour = "black"),
    plot.tag.position = c(0, 1),
    panel.border = ggplot2::element_blank(),
    legend.position = "none",
    plot.margin = ggplot2::margin(3, 4, 2, 2)
  )

pore_colours <- setNames(
  scales::hue_pal()(length(pore_names)),
  pore_names
)


# ============================================================
# 11. DOMINANT-SHAPE PANEL
# ============================================================

make_dominant_shape_panel <- function(result_obj, row_number) {

  overlay_trace <- result_obj$overlay
  median_trace <- result_obj$median_trace
  summary_tbl <- result_obj$summary

  if (
    nrow(overlay_trace) == 0 ||
    nrow(median_trace) == 0 ||
    !identical(summary_tbl$status[1], "OK")
  ) {

    p <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No dominant shape available",
        size = 2.4
      ) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::theme_void(base_size = 8)

    if (row_number == 1) {
      p <- p + ggplot2::labs(tag = "e")
    }

    return(p)
  }

  x_left <- min(median_trace$time_relative_ms, na.rm = TRUE)
  x_right <- max(median_trace$time_relative_ms, na.rm = TRUE)

  y_min <- min(c(overlay_trace$current_centered_pA, median_trace$current_centered_pA, 0), na.rm = TRUE)
  y_max <- max(c(overlay_trace$current_centered_pA, median_trace$current_centered_pA, 0), na.rm = TRUE)

  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range <= 0) {
    y_range <- DOM_CURRENT_SCALEBAR_PA
  }

  y_bottom <- y_min - 0.18 * y_range
  y_top <- max(y_max + 0.20 * y_range, DOM_CURRENT_SCALEBAR_PA * 0.2)

  x_bar_right <- x_right - 0.07 * (x_right - x_left)
  x_bar_left <- x_bar_right - DOM_TIME_SCALEBAR_MS
  y_bar_bottom <- y_bottom + 0.10 * (y_top - y_bottom)
  y_bar_top <- y_bar_bottom + DOM_CURRENT_SCALEBAR_PA

  label_text <- paste0("n = ", summary_tbl$n_in_dominant_cluster[1])

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "grey85",
      linewidth = 0.35
    ) +
    ggplot2::geom_line(
      data = overlay_trace,
      ggplot2::aes(
        x = time_relative_ms,
        y = current_centered_pA,
        group = event_id
      ),
      colour = "grey75",
      linewidth = 0.35,
      alpha = 0.75,
      lineend = "round",
      linejoin = "round"
    ) +
    ggplot2::geom_line(
      data = median_trace,
      ggplot2::aes(
        x = time_relative_ms,
        y = current_centered_pA
      ),
      colour = "black",
      linewidth = 0.70,
      lineend = "round",
      linejoin = "round"
    ) +
    ggplot2::annotate(
      "text",
      x = x_left + 0.10 * (x_right - x_left),
      y = y_top - 0.06 * (y_top - y_bottom),
      label = label_text,
      hjust = 0,
      size = 2.3
    ) +
    ggplot2::annotate(
      "segment",
      x = x_bar_left,
      xend = x_bar_right,
      y = y_bar_bottom,
      yend = y_bar_bottom,
      linewidth = 0.42,
      colour = "black"
    ) +
    ggplot2::annotate(
      "segment",
      x = x_bar_right,
      xend = x_bar_right,
      y = y_bar_bottom,
      yend = y_bar_top,
      linewidth = 0.42,
      colour = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = (x_bar_left + x_bar_right) / 2,
      y = y_bar_bottom - 0.03 * (y_top - y_bottom),
      label = paste0(sprintf("%.1f", DOM_TIME_SCALEBAR_MS), " ms"),
      vjust = 1,
      size = 2.1
    ) +
    ggplot2::annotate(
      "text",
      x = x_bar_right + 0.035 * (x_right - x_left),
      y = (y_bar_bottom + y_bar_top) / 2,
      label = paste0(sprintf("%.0f", DOM_CURRENT_SCALEBAR_PA), " pA"),
      angle = 90,
      hjust = 0.5,
      vjust = 0,
      size = 2.1
    ) +
    ggplot2::coord_cartesian(
      xlim = c(x_left, x_right),
      ylim = c(y_bottom, y_top),
      clip = "off"
    ) +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(2, 5, 2, 6),
      plot.tag = ggplot2::element_text(size = 11, face = "bold", colour = "black"),
      plot.tag.position = c(0, 1)
    )

  if (row_number == 1) {
    p <- p + ggplot2::labs(tag = "e")
  }

  p
}


# ============================================================
# 12. BUILD ONE ROW PER PORE
# ============================================================

make_pore_row <- function(pore_name, row_number) {

  pore_data <- events_plot |>
    dplyr::filter(pore_id == pore_name)

  dwell_hist_data <- pore_data |>
    dplyr::filter(
      dwell_total_ms >= MIN_DWELL_MS,
      dwell_total_ms <= DISPLAY_MAX_DWELL_HIST
    )

  scatter_data <- if (
    is.finite(SCATTER_MAX_POINTS) &&
    nrow(pore_data) > SCATTER_MAX_POINTS
  ) {
    pore_data |>
      dplyr::slice_sample(n = SCATTER_MAX_POINTS)
  } else {
    pore_data
  }

  trace_data <- trace_all |>
    dplyr::filter(pore_id == pore_name)

  dominant_result <- dominant_shape_results[[pore_name]]

  n_events_plotted <- nrow(pore_data)
  this_colour <- unname(pore_colours[pore_name])

  dwell_curve <- fit_shifted_exponential(
    x_ms = dwell_hist_data$dwell_total_ms,
    threshold_ms = MIN_DWELL_MS
  )

  current_curve <- scale_density_to_counts(
    x = pore_data$I_over_I0_pct,
    binwidth = IIO_BINWIDTH,
    x_from = MIN_I_OVER_I0,
    x_to = MAX_I_OVER_I0,
    n_events = nrow(pore_data)
  )

  trace_start_s <- min(trace_data$time_s, na.rm = TRUE)
  trace_end_s <- max(trace_data$time_s, na.rm = TRUE)
  trace_breaks <- pretty(c(trace_start_s, trace_end_s), n = 4)

  baseline_for_plot <- stats::median(trace_data$baseline_pA, na.rm = TRUE)
  trace_y_min <- min(trace_data$current_pA, baseline_for_plot, na.rm = TRUE)
  trace_y_max <- max(trace_data$current_pA, baseline_for_plot, na.rm = TRUE)
  trace_y_range <- trace_y_max - trace_y_min
  if (!is.finite(trace_y_range) || trace_y_range <= 0) {
    trace_y_range <- 10
  }
  trace_limits <- c(trace_y_min - 0.08 * trace_y_range, trace_y_max + 0.08 * trace_y_range)

  scatter_plot <- ggplot2::ggplot(
    scatter_data,
    ggplot2::aes(x = I_over_I0_pct, y = dwell_total_ms)
  ) +
    ggplot2::geom_point(
      colour = this_colour,
      alpha = 0.50,
      size = 0.75
    ) +
    ggplot2::scale_x_continuous(
      limits = c(MIN_I_OVER_I0, MAX_I_OVER_I0),
      breaks = seq(0, 100, by = 20),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_log10(
      limits = c(MIN_DWELL_MS, MAX_DWELL_MS),
      breaks = c(0.2, 0.5, 1, 2, 5, 10, 20, 50, 100),
      labels = c("0.2", "0.5", "1", "2", "5", "10", "20", "50", "100"),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = paste0(pore_name, "   n = ", n_events_plotted),
      x = expression(I/I[0] ~ " (%)"),
      y = "Dwell time (ms)"
    ) +
    article_theme

  dwell_histogram <- ggplot2::ggplot(
    dwell_hist_data,
    ggplot2::aes(x = dwell_total_ms)
  ) +
    ggplot2::geom_histogram(
      binwidth = DWELL_BINWIDTH,
      boundary = MIN_DWELL_MS,
      fill = this_colour,
      colour = this_colour,
      alpha = 0.55,
      linewidth = 0.07
    ) +
    ggplot2::geom_line(
      data = dwell_curve,
      ggplot2::aes(x = x, y = scaled_count),
      inherit.aes = FALSE,
      colour = "grey10",
      linewidth = 0.42
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0.2, 0.5, 1, 2, 3, 4, 5),
      labels = c("0.2", "0.5", "1", "2", "3", "4", "5"),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    ggplot2::coord_cartesian(
      xlim = c(MIN_DWELL_MS, DISPLAY_MAX_DWELL_HIST),
      ylim = c(0, dwell_hist_ymax)
    ) +
    ggplot2::labs(
      x = "Dwell time (ms)",
      y = "Count (N)"
    ) +
    article_theme

  current_histogram <- ggplot2::ggplot(
    pore_data,
    ggplot2::aes(x = I_over_I0_pct)
  ) +
    ggplot2::geom_histogram(
      binwidth = IIO_BINWIDTH,
      boundary = MIN_I_OVER_I0,
      fill = this_colour,
      colour = this_colour,
      alpha = 0.55,
      linewidth = 0.07
    ) +
    ggplot2::geom_line(
      data = current_curve,
      ggplot2::aes(x = x, y = scaled_count),
      inherit.aes = FALSE,
      colour = "grey10",
      linewidth = 0.45
    ) +
    ggplot2::scale_x_continuous(
      limits = c(MIN_I_OVER_I0, MAX_I_OVER_I0),
      breaks = seq(0, 100, by = 20),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, current_hist_ymax)) +
    ggplot2::labs(
      x = expression(I/I[0] ~ " (%)"),
      y = "Count (N)"
    ) +
    article_theme

  trace_plot <- ggplot2::ggplot(
    trace_data,
    ggplot2::aes(x = time_s, y = current_pA)
  ) +
    ggplot2::geom_hline(
      yintercept = baseline_for_plot,
      colour = "grey40",
      linewidth = 0.30,
      linetype = "dashed"
    ) +
    ggplot2::geom_line(
      colour = this_colour,
      linewidth = 0.10
    ) +
    ggplot2::scale_x_continuous(
      limits = c(trace_start_s, trace_end_s),
      breaks = trace_breaks[trace_breaks >= trace_start_s & trace_breaks <= trace_end_s],
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = trace_limits,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      x = "Time (s)",
      y = "Current (pA)"
    ) +
    article_theme

  dominant_shape_plot <- make_dominant_shape_panel(
    result_obj = dominant_result,
    row_number = row_number
  )

  if (row_number == 1) {
    scatter_plot <- scatter_plot + ggplot2::labs(tag = "a")
    dwell_histogram <- dwell_histogram + ggplot2::labs(tag = "b")
    current_histogram <- current_histogram + ggplot2::labs(tag = "c")
    trace_plot <- trace_plot + ggplot2::labs(tag = "d")
  }

  patchwork::wrap_plots(
    scatter_plot,
    dwell_histogram,
    current_histogram,
    trace_plot,
    dominant_shape_plot,
    ncol = 5,
    widths = c(1.12, 1, 1, 1.10, 0.90)
  )
}


# ============================================================
# 13. FINAL FIGURE
# ============================================================

message("Création de la figure...")

pore_rows <- purrr::map2(
  pore_names,
  seq_along(pore_names),
  make_pore_row
)

final_plot <- patchwork::wrap_plots(
  pore_rows,
  ncol = 1
) +
  patchwork::plot_annotation(
    title = paste0(peptide_id, " (", n_pores, " pores)"),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
    )
  )

png_file <- file.path(
  output_folder,
  paste0(peptide_id, "_dominant_shape_only.png")
)

figure_height <- max(
  4,
  HEIGHT_PER_PORE_IN * n_pores + 0.7
)

ggplot2::ggsave(
  filename = png_file,
  plot = final_plot,
  width = FIGURE_WIDTH_IN,
  height = figure_height,
  units = "in",
  dpi = PNG_DPI,
  bg = "white",
  limitsize = FALSE
)

message("Terminé.")
message("PNG créé : ", png_file)
message("Dossier de sortie : ", output_folder)
