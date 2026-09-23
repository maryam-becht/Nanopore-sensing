
required_packages <- c(
  "openxlsx",
  "dplyr",
  "purrr",
  "ggplot2",
  "patchwork",
  "tibble"
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

peptide_id <- basename(
  peptide_folder
)

output_folder <- file.path(
  peptide_folder,
  paste0(
    peptide_id,
    "_R_folder_output"
  )
)

if (!dir.exists(output_folder)) {
  dir.create(
    output_folder,
    recursive = TRUE
  )
}

message(
  "Dossier peptide : ",
  peptide_folder
)

message(
  "Peptide : ",
  peptide_id
)


# ============================================================
# 2. PARAMETERS
# ============================================================

MIN_DWELL_MS <- 0.2
MAX_DWELL_MS <- 100

MIN_I_OVER_I0 <- 0
MAX_I_OVER_I0 <- 100

DISPLAY_MAX_DWELL_HIST <- 5

# Remove pores with fewer than 100 plotted events
MIN_EVENTS_PER_PORE <- 100

DWELL_BINWIDTH <- 0.20
IIO_BINWIDTH <- 2

SCATTER_MAX_POINTS <- 5000
TRACE_MAX_POINTS <- 4500
GMM_MAX_POINTS <- 4000

PNG_DPI <- 180
FIGURE_WIDTH_IN <- 12.5
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
  !grepl(
    "^~\\$",
    basename(excel_files)
  )
]

excel_files <- sort(
  excel_files
)

if (length(excel_files) == 0) {
  stop(
    "Aucun fichier Excel trouvé dans le dossier : ",
    peptide_folder
  )
}

message(
  "Nombre de fichiers pore trouvés : ",
  length(excel_files)
)

print(
  basename(excel_files)
)


# ============================================================
# 4. HELPERS
# ============================================================

extract_pore_id_from_filename <- function(file_path) {

  fname <- tools::file_path_sans_ext(
    basename(file_path)
  )

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

    m <- regexpr(
      "P[0-9]+",
      fname,
      ignore.case = TRUE
    )

    if (m[1] > 0) {
      pore_id <- regmatches(
        fname,
        m
      )
    } else {
      pore_id <- fname
    }
  }

  toupper(pore_id)
}


get_header <- function(wb, sheet_name) {

  header_data <- openxlsx::read.xlsx(
    wb,
    sheet = sheet_name,
    rows = 1,
    colNames = FALSE,
    check.names = FALSE
  )

  as.character(
    unlist(
      header_data[
        1,
        ,
        drop = TRUE
      ]
    )
  )
}


get_sheet_last_row <- function(
  wb,
  sheet_name
) {

  first_column <- openxlsx::read.xlsx(
    wb,
    sheet = sheet_name,
    cols = 1,
    colNames = TRUE,
    check.names = FALSE
  )

  nrow(first_column) + 1L
}


fit_one_gaussian <- function(x) {

  mu <- mean(x)
  sigma <- stats::sd(x)

  if (!is.finite(sigma) || sigma < 1e-6) {
    sigma <- 1e-6
  }

  log_likelihood <- sum(
    stats::dnorm(
      x,
      mean = mu,
      sd = sigma,
      log = TRUE
    )
  )

  bic <- -2 * log_likelihood +
    2 * log(length(x))

  list(
    G = 1L,
    means = mu,
    sigmas = sigma,
    proportions = 1,
    bic = bic
  )
}


fit_two_gaussians_em <- function(
  x,
  max_iter = 150,
  tolerance = 1e-6
) {

  n <- length(x)

  means <- as.numeric(
    stats::quantile(
      x,
      probs = c(0.33, 0.67),
      names = FALSE,
      type = 7
    )
  )

  overall_sd <- stats::sd(x)

  if (!is.finite(overall_sd) || overall_sd < 1e-6) {
    overall_sd <- 1
  }

  sigmas <- rep(
    max(overall_sd / 2, 0.5),
    2
  )

  proportions <- c(0.5, 0.5)

  previous_log_likelihood <- -Inf

  for (iteration in seq_len(max_iter)) {

    log_component_1 <- log(
      max(proportions[1], 1e-12)
    ) +
      stats::dnorm(
        x,
        mean = means[1],
        sd = sigmas[1],
        log = TRUE
      )

    log_component_2 <- log(
      max(proportions[2], 1e-12)
    ) +
      stats::dnorm(
        x,
        mean = means[2],
        sd = sigmas[2],
        log = TRUE
      )

    row_max <- pmax(
      log_component_1,
      log_component_2
    )

    denominator <- row_max +
      log(
        exp(log_component_1 - row_max) +
          exp(log_component_2 - row_max)
      )

    responsibility_1 <- exp(
      log_component_1 - denominator
    )

    responsibility_2 <- 1 - responsibility_1

    effective_n_1 <- sum(
      responsibility_1
    )

    effective_n_2 <- sum(
      responsibility_2
    )

    if (
      effective_n_1 < 2 ||
      effective_n_2 < 2
    ) {
      return(NULL)
    }

    proportions <- c(
      effective_n_1 / n,
      effective_n_2 / n
    )

    means <- c(
      sum(responsibility_1 * x) / effective_n_1,
      sum(responsibility_2 * x) / effective_n_2
    )

    sigmas <- sqrt(
      c(
        sum(
          responsibility_1 * (x - means[1])^2
        ) / effective_n_1,

        sum(
          responsibility_2 * (x - means[2])^2
        ) / effective_n_2
      )
    )

    sigmas <- pmax(
      sigmas,
      0.25
    )

    log_likelihood <- sum(
      denominator
    )

    if (
      is.finite(previous_log_likelihood) &&
      abs(log_likelihood - previous_log_likelihood) < tolerance
    ) {
      previous_log_likelihood <- log_likelihood
      break
    }

    previous_log_likelihood <- log_likelihood
  }

  peak_order <- order(
    means
  )

  means <- means[peak_order]
  sigmas <- sigmas[peak_order]
  proportions <- proportions[peak_order]

  bic <- -2 * previous_log_likelihood +
    5 * log(n)

  list(
    G = 2L,
    means = means,
    sigmas = sigmas,
    proportions = proportions,
    bic = bic
  )
}


fit_gmm_fast <- function(
  pore_data,
  pore_name
) {

  x <- pore_data$I_over_I0_pct
  x <- x[is.finite(x)]

  if (length(x) < 5) {
    return(
      list(
        curve = tibble::tibble(
          x = numeric(0),
          scaled_count = numeric(0)
        ),
        summary = tibble::tibble(
          pore_id = pore_name,
          n_events_total = length(x),
          n_events_fit = length(x),
          selected_components = NA_integer_,
          selected_model = NA_character_,
          peak = NA_character_,
          mu_I_over_I0_pct = NA_real_,
          sigma_I_over_I0_pct = NA_real_,
          proportion = NA_real_
        )
      )
    )
  }

  x_fit <- if (length(x) > GMM_MAX_POINTS) {
    sample(
      x,
      size = GMM_MAX_POINTS,
      replace = FALSE
    )
  } else {
    x
  }

  fit_G1 <- fit_one_gaussian(
    x_fit
  )

  fit_G2 <- if (length(x_fit) >= 40) {
    fit_two_gaussians_em(
      x_fit
    )
  } else {
    NULL
  }

  if (
    !is.null(fit_G2) &&
    is.finite(fit_G2$bic) &&
    fit_G2$bic < fit_G1$bic
  ) {
    selected_fit <- fit_G2
    selected_model <- "Two-Gaussian EM"
  } else {
    selected_fit <- fit_G1
    selected_model <- "Single Gaussian"
  }

  selected_G <- selected_fit$G
  means <- selected_fit$means
  sigmas <- selected_fit$sigmas
  proportions <- selected_fit$proportions

  x_grid <- seq(
    MIN_I_OVER_I0,
    MAX_I_OVER_I0,
    length.out = 300
  )

  component_curves <- purrr::map_dfr(
    seq_len(selected_G),
    function(j) {

      density_j <- proportions[j] *
        stats::dnorm(
          x_grid,
          mean = means[j],
          sd = sigmas[j]
        )

      tibble::tibble(
        x = x_grid,
        scaled_count =
          density_j *
          length(x) *
          IIO_BINWIDTH
      )
    }
  )

  total_curve <- component_curves |>
    dplyr::group_by(x) |>
    dplyr::summarise(
      scaled_count = sum(scaled_count),
      .groups = "drop"
    )

  summary_table <- tibble::tibble(
    pore_id = pore_name,
    n_events_total = length(x),
    n_events_fit = length(x_fit),
    selected_components = selected_G,
    selected_model = selected_model,
    peak = paste0(
      "Peak ",
      seq_len(selected_G)
    ),
    mu_I_over_I0_pct = means,
    sigma_I_over_I0_pct = sigmas,
    proportion = proportions
  )

  list(
    curve = total_curve,
    summary = summary_table
  )
}


fit_shifted_exponential <- function(
  pore_data,
  pore_name
) {

  dwell_fit <- pore_data$dwell_total_ms[
    is.finite(pore_data$dwell_total_ms) &
      pore_data$dwell_total_ms > MIN_DWELL_MS &
      pore_data$dwell_total_ms <= DISPLAY_MAX_DWELL_HIST
  ]

  shifted_dwell <- dwell_fit - MIN_DWELL_MS

  if (
    length(shifted_dwell) < 5 ||
    mean(shifted_dwell) <= 0
  ) {
    return(
      list(
        curve = tibble::tibble(
          x = numeric(0),
          scaled_count = numeric(0)
        ),
        summary = tibble::tibble(
          pore_id = pore_name,
          n_events_used = length(shifted_dwell),
          lambda_per_ms = NA_real_,
          tau_above_threshold_ms = NA_real_
        )
      )
    )
  }

  lambda_hat <- 1 / mean(
    shifted_dwell
  )

  x_grid <- seq(
    MIN_DWELL_MS,
    DISPLAY_MAX_DWELL_HIST,
    length.out = 250
  )

  density_fit <- lambda_hat *
    exp(
      -lambda_hat * (x_grid - MIN_DWELL_MS)
    )

  list(
    curve = tibble::tibble(
      x = x_grid,
      scaled_count =
        density_fit *
        length(dwell_fit) *
        DWELL_BINWIDTH
    ),
    summary = tibble::tibble(
      pore_id = pore_name,
      n_events_used = length(dwell_fit),
      lambda_per_ms = lambda_hat,
      tau_above_threshold_ms = 1 / lambda_hat
    )
  )
}


read_events_one_file <- function(
  file_path,
  pore_name
) {

  wb <- openxlsx::loadWorkbook(
    file_path
  )

  sheet_names <- names(wb)

  if (!"events" %in% sheet_names) {
    stop(
      "La feuille 'events' manque dans ",
      basename(file_path)
    )
  }

  header <- get_header(
    wb,
    "events"
  )

  dwell_col <- match(
    "dwell_total_ms",
    header
  )

  iio_col <- match(
    "I_over_I0_pct",
    header
  )

  if (is.na(dwell_col) || is.na(iio_col)) {
    stop(
      "Colonnes dwell_total_ms / I_over_I0_pct manquantes dans ",
      basename(file_path)
    )
  }

  df <- openxlsx::read.xlsx(
    wb,
    sheet = "events",
    cols = unique(c(dwell_col, iio_col)),
    colNames = TRUE,
    check.names = FALSE
  )

  df |>
    dplyr::transmute(
      pore_id = pore_name,
      dwell_total_ms = as.numeric(dwell_total_ms),
      I_over_I0_pct = as.numeric(I_over_I0_pct)
    ) |>
    dplyr::filter(
      is.finite(dwell_total_ms),
      is.finite(I_over_I0_pct)
    )
}


read_sampled_trace <- function(
  file_path,
  pore_name
) {

  wb <- openxlsx::loadWorkbook(
    file_path
  )

  sheet_names <- names(wb)
  trace_sheet <- NA_character_

  if ("pore_summary" %in% sheet_names) {

    ps <- tryCatch(
      openxlsx::read.xlsx(
        wb,
        sheet = "pore_summary",
        check.names = FALSE
      ),
      error = function(e) NULL
    )

    if (
      !is.null(ps) &&
      "trace_sheet" %in% names(ps) &&
      length(ps$trace_sheet) >= 1
    ) {
      ts <- as.character(ps$trace_sheet[1])
      if (!is.na(ts) && ts != "" && ts %in% sheet_names) {
        trace_sheet <- ts
      }
    }
  }

  if (is.na(trace_sheet)) {
    trace_candidates <- sheet_names[
      grepl("trace", sheet_names, ignore.case = TRUE)
    ]
    if (length(trace_candidates) >= 1) {
      trace_sheet <- trace_candidates[1]
    }
  }

  if (is.na(trace_sheet)) {
    stop(
      "Aucune feuille de trace trouvée dans ",
      basename(file_path)
    )
  }

  header <- get_header(
    wb,
    trace_sheet
  )

  time_col <- match("time_s", header)
  current_col <- match("current_pA", header)
  baseline_col <- match("baseline_pA", header)

  if (is.na(time_col) || is.na(current_col)) {
    stop(
      "Colonnes time_s/current_pA manquantes dans ",
      basename(file_path),
      " / ",
      trace_sheet
    )
  }

  selected_columns <- c(time_col, current_col)
  if (!is.na(baseline_col)) {
    selected_columns <- c(selected_columns, baseline_col)
  }

  last_row <- get_sheet_last_row(
    wb,
    trace_sheet
  )

  data_rows <- seq.int(
    from = 2,
    to = last_row
  )

  if (length(data_rows) > TRACE_MAX_POINTS) {
    data_rows <- unique(
      round(
        seq(
          2,
          last_row,
          length.out = TRACE_MAX_POINTS
        )
      )
    )
  }

  selected_rows <- c(1, data_rows)

  trace_data <- openxlsx::read.xlsx(
    wb,
    sheet = trace_sheet,
    rows = selected_rows,
    cols = selected_columns,
    colNames = TRUE,
    check.names = FALSE
  )

  baseline_value <- if (
    "baseline_pA" %in% names(trace_data)
  ) {
    stats::median(
      as.numeric(trace_data$baseline_pA),
      na.rm = TRUE
    )
  } else {
    stats::median(
      as.numeric(trace_data$current_pA),
      na.rm = TRUE
    )
  }

  trace_data |>
    dplyr::transmute(
      pore_id = pore_name,
      time_s = as.numeric(time_s),
      current_pA = as.numeric(current_pA),
      baseline_pA = baseline_value,
      current_iio_pct = (as.numeric(current_pA) / baseline_value) * 100
    ) |>
    dplyr::filter(
      is.finite(time_s),
      is.finite(current_pA),
      is.finite(current_iio_pct)
    )
}


# ============================================================
# 5. READ ALL PORE FILES
# ============================================================

pore_ids <- vapply(
  excel_files,
  extract_pore_id_from_filename,
  character(1)
)

pore_info <- tibble::tibble(
  file_path = excel_files,
  file_name = basename(excel_files),
  pore_id = pore_ids
) |>
  dplyr::mutate(
    pore_number = suppressWarnings(
      as.integer(
        sub("^P", "", pore_id)
      )
    )
  ) |>
  dplyr::arrange(
    pore_number,
    pore_id
  ) |>
  dplyr::distinct(
    pore_id,
    .keep_all = TRUE
  )

pore_names <- pore_info$pore_id
n_pores <- nrow(pore_info)

message(
  "Nombre de pores détectés : ",
  n_pores
)

print(
  pore_info |>
    dplyr::select(pore_id, file_name)
)

events_raw <- purrr::map2_dfr(
  pore_info$file_path,
  pore_info$pore_id,
  read_events_one_file
)

events_plot <- events_raw |>
  dplyr::filter(
    dwell_total_ms >= MIN_DWELL_MS,
    dwell_total_ms <= MAX_DWELL_MS,
    I_over_I0_pct >= MIN_I_OVER_I0,
    I_over_I0_pct <= MAX_I_OVER_I0
  ) |>
  dplyr::mutate(
    pore_id = factor(
      pore_id,
      levels = pore_names
    )
  )

event_counts <- tibble::tibble(
  pore_id = pore_names
) |>
  dplyr::left_join(
    events_raw |>
      dplyr::group_by(pore_id) |>
      dplyr::summarise(
        n_events_in_excel = dplyr::n(),
        .groups = "drop"
      ),
    by = "pore_id"
  ) |>
  dplyr::left_join(
    events_plot |>
      dplyr::mutate(
        pore_id = as.character(pore_id)
      ) |>
      dplyr::group_by(pore_id) |>
      dplyr::summarise(
        n_events_plotted = dplyr::n(),
        .groups = "drop"
      ),
    by = "pore_id"
  ) |>
  dplyr::mutate(
    n_events_in_excel = dplyr::coalesce(n_events_in_excel, 0L),
    n_events_plotted = dplyr::coalesce(n_events_plotted, 0L)
  )

print(event_counts)

# ============================================================
# 5B. REMOVE PORES WITH FEWER THAN 100 PLOTTED EVENTS
# ============================================================

kept_pores <- event_counts |>
  dplyr::filter(
    n_events_plotted >= MIN_EVENTS_PER_PORE
  ) |>
  dplyr::pull(
    pore_id
  )

if (length(kept_pores) == 0) {
  stop(
    "Aucun pore ne passe le filtre de ",
    MIN_EVENTS_PER_PORE,
    " événements."
  )
}

pore_info <- pore_info |>
  dplyr::filter(
    pore_id %in% kept_pores
  )

pore_names <- pore_info$pore_id
n_pores <- length(pore_names)

events_raw <- events_raw |>
  dplyr::filter(
    pore_id %in% kept_pores
  )

events_plot <- events_plot |>
  dplyr::filter(
    as.character(pore_id) %in% kept_pores
  ) |>
  dplyr::mutate(
    pore_id = factor(
      as.character(pore_id),
      levels = pore_names
    )
  )

event_counts <- event_counts |>
  dplyr::filter(
    pore_id %in% kept_pores
  )

message(
  "Pores gardés (n >= ",
  MIN_EVENTS_PER_PORE,
  ") : ",
  paste(kept_pores, collapse = ", ")
)

write.csv(
  event_counts,
  file.path(
    output_folder,
    paste0(peptide_id, "_event_counts.csv")
  ),
  row.names = FALSE
)

write.csv(
  pore_info |>
    dplyr::select(pore_id, file_name),
  file.path(
    output_folder,
    paste0(peptide_id, "_pore_mapping.csv")
  ),
  row.names = FALSE
)

trace_all <- purrr::map2_dfr(
  pore_info$file_path,
  pore_info$pore_id,
  read_sampled_trace
) |>
  dplyr::mutate(
    pore_id = factor(
      pore_id,
      levels = pore_names
    )
  )


# ============================================================
# 6. FITS
# ============================================================

gmm_results <- setNames(
  purrr::map(
    pore_names,
    function(current_pore) {
      fit_gmm_fast(
        pore_data = events_plot |>
          dplyr::filter(pore_id == current_pore),
        pore_name = current_pore
      )
    }
  ),
  pore_names
)

gmm_summary <- purrr::map_dfr(
  gmm_results,
  "summary"
)

write.csv(
  gmm_summary,
  file.path(
    output_folder,
    paste0(peptide_id, "_Gaussian_summary.csv")
  ),
  row.names = FALSE
)

exponential_results <- setNames(
  purrr::map(
    pore_names,
    function(current_pore) {
      fit_shifted_exponential(
        pore_data = events_plot |>
          dplyr::filter(pore_id == current_pore),
        pore_name = current_pore
      )
    }
  ),
  pore_names
)

exponential_summary <- purrr::map_dfr(
  exponential_results,
  "summary"
)

write.csv(
  exponential_summary,
  file.path(
    output_folder,
    paste0(peptide_id, "_exponential_summary.csv")
  ),
  row.names = FALSE
)


# ============================================================
# 7. COMMON AXES
#
# The two middle panels use Count (N).
# Their y-axis limits are calculated globally and then applied
# identically to every pore.
# ============================================================

dwell_hist_ymax <- max(
  purrr::map_dbl(
    pore_names,
    function(current_pore) {

      pore_data <- events_plot |>
        dplyr::filter(
          pore_id == current_pore,
          dwell_total_ms >= MIN_DWELL_MS,
          dwell_total_ms <= DISPLAY_MAX_DWELL_HIST
        )

      hist_max <- if (nrow(pore_data) > 0) {

        dwell_breaks_for_count <- seq(
          MIN_DWELL_MS,
          DISPLAY_MAX_DWELL_HIST + DWELL_BINWIDTH,
          by = DWELL_BINWIDTH
        )

        h <- hist(
          pore_data$dwell_total_ms,
          breaks = dwell_breaks_for_count,
          plot = FALSE,
          include.lowest = TRUE,
          right = FALSE
        )

        max(
          h$counts,
          na.rm = TRUE
        )

      } else {
        0
      }

      curve_max <- if (
        nrow(
          exponential_results[[current_pore]]$curve
        ) > 0
      ) {

        max(
          exponential_results[[current_pore]]$curve$scaled_count,
          na.rm = TRUE
        )

      } else {
        0
      }

      max(
        hist_max,
        curve_max
      )
    }
  ),
  na.rm = TRUE
) * 1.05


current_hist_ymax <- max(
  purrr::map_dbl(
    pore_names,
    function(current_pore) {

      pore_data <- events_plot |>
        dplyr::filter(
          pore_id == current_pore
        )

      hist_max <- if (nrow(pore_data) > 0) {

        current_breaks_for_count <- seq(
          MIN_I_OVER_I0,
          MAX_I_OVER_I0 + IIO_BINWIDTH,
          by = IIO_BINWIDTH
        )

        h <- hist(
          pore_data$I_over_I0_pct,
          breaks = current_breaks_for_count,
          plot = FALSE,
          include.lowest = TRUE,
          right = FALSE
        )

        max(
          h$counts,
          na.rm = TRUE
        )

      } else {
        0
      }

      curve_max <- if (
        nrow(
          gmm_results[[current_pore]]$curve
        ) > 0
      ) {

        max(
          gmm_results[[current_pore]]$curve$scaled_count,
          na.rm = TRUE
        )

      } else {
        0
      }

      max(
        hist_max,
        curve_max
      )
    }
  ),
  na.rm = TRUE
) * 1.05


# Same RAW current scale (pA) for all retained pores.
# Each pore therefore keeps its real baseline level instead of being forced to 100%.
trace_min <- min(
  trace_all$current_pA,
  na.rm = TRUE
)

trace_max <- max(
  trace_all$current_pA,
  na.rm = TRUE
)

trace_range <- trace_max - trace_min

if (!is.finite(trace_range) || trace_range <= 0) {
  trace_range <- max(
    abs(trace_max),
    1
  )
}

trace_margin <- 0.05 * trace_range

TRACE_Y_LIMITS <- c(
  trace_min - trace_margin,
  trace_max + trace_margin
)


# ============================================================
# 8. COLOURS + THEME
# ============================================================

colorblind_palette <- c(
  "#4E79A7",  # muted blue
  "#E07B62",  # soft coral
  "#6F9A78",  # softer muted green
  "#9C6FAE",  # soft purple
  "#6B6B6B",  # neutral grey
  "#76A5AF",  # blue-grey
  "#B07A5A",  # muted terracotta
  "#7E6E8E"   # grey-purple
)

pore_colours <- rep(
  colorblind_palette,
  length.out = n_pores
)

names(pore_colours) <- pore_names

article_theme <- ggplot2::theme_classic(base_size = 8) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(colour = "black", size = 6.7),
    axis.title = ggplot2::element_text(colour = "black", size = 7.7),
    plot.title = ggplot2::element_text(colour = "black", size = 8.7, hjust = 0.5),
    plot.tag = ggplot2::element_text(colour = "black", size = 10.5, face = "bold"),
    plot.tag.position = c(0, 1),
    axis.line = ggplot2::element_line(colour = "black", linewidth = 0.30),
    axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.25),
    legend.position = "none",
    plot.margin = ggplot2::margin(2, 3, 2, 3)
  )


# ============================================================
# 9. BUILD ONE ROW PER PORE
# ============================================================

make_pore_row <- function(
  pore_name,
  row_number
) {

  pore_data <- events_plot |>
    dplyr::filter(pore_id == pore_name)

  # Panel b uses only events displayed between 0.2 and 5 ms
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
      dplyr::slice_sample(
        n = SCATTER_MAX_POINTS
      )
  } else {
    pore_data
  }

  trace_data <- trace_all |>
    dplyr::filter(pore_id == pore_name)

  n_events_plotted <- nrow(pore_data)
  this_colour <- unname(pore_colours[pore_name])
  gmm_curve <- gmm_results[[pore_name]]$curve
  exponential_curve <- exponential_results[[pore_name]]$curve

  trace_start_s <- min(trace_data$time_s, na.rm = TRUE)
  trace_end_s <- max(trace_data$time_s, na.rm = TRUE)

  trace_breaks <- pretty(
    c(trace_start_s, trace_end_s),
    n = 4
  )

  scatter_plot <- ggplot2::ggplot(
    scatter_data,
    ggplot2::aes(
      x = I_over_I0_pct,
      y = dwell_total_ms
    )
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
      x = expression(I/I[0]~" (%)"),
      y = "Dwell time (ms)"
    ) +
    article_theme

  dwell_histogram <- ggplot2::ggplot(
    dwell_hist_data,
    ggplot2::aes(
      x = dwell_total_ms
    )
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
      data = exponential_curve,
      ggplot2::aes(
        x = x,
        y = scaled_count
      ),
      inherit.aes = FALSE,
      colour = "grey10",
      linewidth = 0.42
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(
        0.2,
        0.5,
        1,
        2,
        3,
        4,
        5
      ),
      labels = c(
        "0.2",
        "0.5",
        "1",
        "2",
        "3",
        "4",
        "5"
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.03
        )
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        MIN_DWELL_MS,
        DISPLAY_MAX_DWELL_HIST
      ),
      ylim = c(
        0,
        dwell_hist_ymax
      )
    ) +
    ggplot2::labs(
      x = "Dwell time (ms)",
      y = "Count (N)"
    ) +
    article_theme

  current_histogram <- ggplot2::ggplot(
    pore_data,
    ggplot2::aes(
      x = I_over_I0_pct
    )
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
      data = gmm_curve,
      ggplot2::aes(
        x = x,
        y = scaled_count
      ),
      inherit.aes = FALSE,
      colour = "grey10",
      linewidth = 0.45
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        MIN_I_OVER_I0,
        MAX_I_OVER_I0
      ),
      breaks = seq(
        0,
        100,
        by = 20
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.03
        )
      )
    ) +
    ggplot2::coord_cartesian(
      ylim = c(
        0,
        current_hist_ymax
      )
    ) +
    ggplot2::labs(
      x = expression(
        I/I[0]~" (%)"
      ),
      y = "Count (N)"
    ) +
    article_theme

  # Real baseline for this pore, in pA
  baseline_for_plot <- stats::median(
    trace_data$baseline_pA,
    na.rm = TRUE
  )

  trace_plot <- ggplot2::ggplot(
    trace_data,
    ggplot2::aes(
      x = time_s,
      y = current_pA
    )
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
      limits = c(
        trace_start_s,
        trace_end_s
      ),
      breaks = trace_breaks[
        trace_breaks >= trace_start_s &
          trace_breaks <= trace_end_s
      ],
      expand = ggplot2::expansion(
        mult = c(0, 0)
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = TRACE_Y_LIMITS,
      expand = ggplot2::expansion(
        mult = c(0, 0)
      )
    ) +
    ggplot2::labs(
      x = "Time (s)",
      y = "Current (pA)"
    ) +
    article_theme

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
    ncol = 4,
    widths = c(1.12, 1, 1, 1.12)
  )
}


# ============================================================
# 10. FINAL FIGURE
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
  paste0(
    peptide_id,
    "_counts_dwell5ms_raw_current_true_baseline.png"
  )
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
