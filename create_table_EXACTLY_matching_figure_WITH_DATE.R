# ============================================================
# TABLE QUI CORRESPOND EXACTEMENT AUX FIGURES
# ============================================================
#
# IMPORTANT :
# Cette table reproduit exactement la logique du script de figure :
#
#   dwell_total_ms >= 0.2
#   dwell_total_ms <= 100
#   I_over_I0_pct >= 0
#   I_over_I0_pct <= 100
#   garder les pores avec n_events_plotted >= 100
#
# Le numero/label du pore correspond EXACTEMENT au P1, P2, P3...
# affiche sur la figure, car il est extrait du nom du fichier comme
# dans le script de figure.
#
# Un dossier peptide a la fois, exactement comme le script de figure.
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

required_packages <- c(
  "openxlsx",
  "dplyr",
  "purrr",
  "tibble",
  "stringr"
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
    caption = "Choisis le dossier du peptide utilise pour la figure"
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
  stop("Aucun dossier selectionne.")
}

peptide_folder <- normalizePath(
  peptide_folder,
  winslash = "/",
  mustWork = TRUE
)

peptide_id <- basename(peptide_folder)

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


# ============================================================
# 2. EXACT FIGURE PARAMETERS
# ============================================================

MIN_DWELL_MS <- 0.2
MAX_DWELL_MS <- 100

MIN_I_OVER_I0 <- 0
MAX_I_OVER_I0 <- 100

# EXACTEMENT comme le code de figure :
# un pore est garde si n_events_plotted >= 100.
MIN_EVENTS_PER_PORE <- 100

TRACE_MAX_POINTS <- 4500


# ============================================================
# 3. FIND EXCEL FILES EXACTLY LIKE FIGURE SCRIPT
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

# Ne pas lire les tables recap creees par ce script.
excel_files <- excel_files[
  !grepl(
    "_figure_matching_table\\.xlsx$",
    basename(excel_files),
    ignore.case = TRUE
  )
]

excel_files <- sort(excel_files)

if (length(excel_files) == 0) {
  stop(
    "Aucun fichier Excel trouve dans : ",
    peptide_folder
  )
}


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


as_numeric_safe <- function(x) {

  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  x <- trimws(as.character(x))
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("%", "", x, fixed = TRUE)

  suppressWarnings(
    as.numeric(x)
  )
}


first_finite <- function(x) {

  x <- as_numeric_safe(x)
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    NA_real_
  } else {
    x[1]
  }
}


parse_date_safe <- function(x) {

  if (length(x) == 0 || all(is.na(x))) {
    return(as.Date(NA))
  }

  x <- x[which(!is.na(x))[1]]

  if (inherits(x, "Date")) {
    return(as.Date(x))
  }

  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }

  if (is.numeric(x) && is.finite(x)) {
    return(
      tryCatch(
        as.Date(
          x,
          origin = "1899-12-30"
        ),
        error = function(e) as.Date(NA)
      )
    )
  }

  txt <- trimws(as.character(x))

  formats <- c(
    "%Y-%m-%d",
    "%d-%m-%Y",
    "%d/%m/%Y",
    "%d-%m-%y",
    "%d/%m/%y",
    "%Y_%m_%d",
    "%d_%m_%Y",
    "%d_%m_%y"
  )

  for (fmt in formats) {

    ans <- tryCatch(
      suppressWarnings(
        as.Date(
          txt,
          format = fmt
        )
      ),
      error = function(e) as.Date(NA)
    )

    if (!is.na(ans)) {
      return(ans)
    }
  }

  as.Date(NA)
}


extract_date_from_text <- function(text) {

  text <- paste(
    as.character(text),
    collapse = " "
  )

  # dd-mm-yy / dd_mm_yy / dd-mm-yyyy
  m <- stringr::str_match(
    text,
    "(?<![0-9])([0-3]?[0-9])[-_/]([01]?[0-9])[-_/]([0-9]{4}|[0-9]{2})(?![0-9])"
  )

  if (!is.na(m[1, 1])) {

    day <- as.integer(m[1, 2])
    month <- as.integer(m[1, 3])
    year <- as.integer(m[1, 4])

    if (year < 100) {
      year <- 2000 + year
    }

    return(
      suppressWarnings(
        as.Date(
          sprintf(
            "%04d-%02d-%02d",
            year,
            month,
            day
          )
        )
      )
    )
  }

  as.Date(NA)
}


# ============================================================
# 5. READ EXACT PLOTTED EVENTS
# ============================================================

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

  if (
    is.na(dwell_col) ||
    is.na(iio_col)
  ) {
    stop(
      "Colonnes dwell_total_ms / I_over_I0_pct manquantes dans ",
      basename(file_path)
    )
  }

  df <- openxlsx::read.xlsx(
    wb,
    sheet = "events",
    cols = unique(
      c(
        dwell_col,
        iio_col
      )
    ),
    colNames = TRUE,
    check.names = FALSE
  )

  df |>
    dplyr::transmute(
      pore_id = pore_name,
      dwell_total_ms =
        as.numeric(dwell_total_ms),
      I_over_I0_pct =
        as.numeric(I_over_I0_pct)
    ) |>
    dplyr::filter(
      is.finite(dwell_total_ms),
      is.finite(I_over_I0_pct)
    )
}


# ============================================================
# 6. READ TRACE EXACTLY LIKE FIGURE + KEEP BASELINE RANGE
# ============================================================

read_trace_summary <- function(
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

      ts <- as.character(
        ps$trace_sheet[1]
      )

      if (
        !is.na(ts) &&
        ts != "" &&
        ts %in% sheet_names
      ) {
        trace_sheet <- ts
      }
    }

  } else {
    ps <- NULL
  }

  if (is.na(trace_sheet)) {

    trace_candidates <- sheet_names[
      grepl(
        "trace",
        sheet_names,
        ignore.case = TRUE
      )
    ]

    if (length(trace_candidates) >= 1) {
      trace_sheet <- trace_candidates[1]
    }
  }

  if (is.na(trace_sheet)) {
    stop(
      "Aucune feuille de trace trouvee dans ",
      basename(file_path)
    )
  }

  header <- get_header(
    wb,
    trace_sheet
  )

  time_col <- match(
    "time_s",
    header
  )

  current_col <- match(
    "current_pA",
    header
  )

  baseline_col <- match(
    "baseline_pA",
    header
  )

  if (
    is.na(time_col) ||
    is.na(current_col)
  ) {
    stop(
      "Colonnes time_s/current_pA manquantes dans ",
      basename(file_path),
      " / ",
      trace_sheet
    )
  }

  selected_columns <- c(
    time_col,
    current_col
  )

  if (!is.na(baseline_col)) {
    selected_columns <- c(
      selected_columns,
      baseline_col
    )
  }

  last_row <- get_sheet_last_row(
    wb,
    trace_sheet
  )

  data_rows <- seq.int(
    from = 2,
    to = last_row
  )

  # EXACTEMENT le meme sous-echantillonnage que la figure.
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

  selected_rows <- c(
    1,
    data_rows
  )

  trace_data <- openxlsx::read.xlsx(
    wb,
    sheet = trace_sheet,
    rows = selected_rows,
    cols = selected_columns,
    colNames = TRUE,
    check.names = FALSE
  )

  current_values <- as_numeric_safe(
    trace_data$current_pA
  )

  if (
    "baseline_pA" %in%
      names(trace_data)
  ) {

    baseline_values <- as_numeric_safe(
      trace_data$baseline_pA
    )

    baseline_values <- baseline_values[
      is.finite(baseline_values)
    ]

  } else {

    baseline_values <- numeric()
  }

  # EXACT baseline_for_plot :
  # le script de figure utilise median(baseline_pA)
  # s'il existe, sinon median(current_pA).
  if (length(baseline_values) > 0) {

    baseline_for_plot <- stats::median(
      baseline_values,
      na.rm = TRUE
    )

    baseline_min <- min(
      baseline_values,
      na.rm = TRUE
    )

    baseline_max <- max(
      baseline_values,
      na.rm = TRUE
    )

  } else {

    finite_current <- current_values[
      is.finite(current_values)
    ]

    baseline_for_plot <- stats::median(
      finite_current,
      na.rm = TRUE
    )

    # Si aucune vraie colonne baseline_pA n'existe,
    # il n'est pas correct d'appeler tout le courant une
    # "baseline range". On laisse donc la range vide.
    baseline_min <- NA_real_
    baseline_max <- NA_real_
  }

  # Duration : prefer pore_summary analyzed duration.
  duration_s <- NA_real_

  if (!is.null(ps)) {

    duration_candidates <- c(
      "analyzed_duration_s",
      "total_analyzed_duration_s",
      "duration_s",
      "recording_duration_s",
      "pore_duration_s"
    )

    hit <- duration_candidates[
      duration_candidates %in% names(ps)
    ]

    if (length(hit) > 0) {
      duration_s <- first_finite(
        ps[[hit[1]]]
      )
    }

    if (!is.finite(duration_s)) {

      duration_min_candidates <- c(
        "analyzed_duration_min",
        "duration_min",
        "recording_duration_min",
        "pore_duration_min"
      )

      hit_min <- duration_min_candidates[
        duration_min_candidates %in% names(ps)
      ]

      if (length(hit_min) > 0) {
        duration_s <- 60 *
          first_finite(
            ps[[hit_min[1]]]
          )
      }
    }
  }

  # Fallback from trace time span.
  if (!is.finite(duration_s)) {

    finite_time <- as_numeric_safe(
      trace_data$time_s
    )

    finite_time <- finite_time[
      is.finite(finite_time)
    ]

    if (length(finite_time) >= 2) {
      duration_s <-
        max(finite_time) -
        min(finite_time)
    }
  }

  # Lipid concentration and date from pore_summary if present.
  lipid_concentration <- NA_real_
  pore_date <- as.Date(NA)

  if (!is.null(ps)) {

    lipid_candidates <- c(
      "lipid_concentration_mg_mL",
      "lipid_concentration_mg_ml",
      "lipid_concentration",
      "DPhPC_concentration_mg_mL",
      "dphpc_concentration_mg_ml"
    )

    lipid_hit <- lipid_candidates[
      lipid_candidates %in% names(ps)
    ]

    if (length(lipid_hit) > 0) {
      lipid_concentration <-
        first_finite(
          ps[[lipid_hit[1]]]
        )
    }

    # 1) Chercher d'abord une vraie colonne date dans pore_summary.
    date_candidates <- c(
      "date",
      "Date",
      "recording_date",
      "Recording_date",
      "experiment_date",
      "Experiment_date"
    )

    date_hit <- date_candidates[
      date_candidates %in% names(ps)
    ]

    if (length(date_hit) > 0) {
      pore_date <- parse_date_safe(
        ps[[date_hit[1]]]
      )
    }

    # 2) Si pas de colonne date, chercher la date dans les champs texte
    #    du pore_summary, notamment source_pore_folder.
    if (is.na(pore_date)) {

      text_candidates <- c(
        "source_pore_folder",
        "source_folder",
        "recording_folder",
        "folder",
        "relative_path",
        "source_file",
        "file_name"
      )

      text_hit <- text_candidates[
        text_candidates %in% names(ps)
      ]

      if (length(text_hit) > 0) {

        for (current_col in text_hit) {

          candidate_text <- paste(
            as.character(ps[[current_col]]),
            collapse = " "
          )

          candidate_date <- extract_date_from_text(
            candidate_text
          )

          if (!is.na(candidate_date)) {
            pore_date <- candidate_date
            break
          }
        }
      }
    }
  }

  # 3) Chercher ensuite dans les autres feuilles du classeur
  #    (file_summary, parameters, README, etc.).
  if (is.na(pore_date)) {

    for (candidate_sheet in sheet_names) {

      candidate_data <- tryCatch(
        openxlsx::read.xlsx(
          wb,
          sheet = candidate_sheet,
          check.names = FALSE
        ),
        error = function(e) NULL
      )

      if (
        is.null(candidate_data) ||
        nrow(candidate_data) == 0
      ) {
        next
      }

      # Chercher une colonne date explicite.
      candidate_date_cols <- intersect(
        c(
          "date",
          "Date",
          "recording_date",
          "Recording_date",
          "experiment_date",
          "Experiment_date"
        ),
        names(candidate_data)
      )

      if (length(candidate_date_cols) > 0) {

        for (current_col in candidate_date_cols) {

          candidate_date <- parse_date_safe(
            candidate_data[[current_col]]
          )

          if (!is.na(candidate_date)) {
            pore_date <- candidate_date
            break
          }
        }
      }

      if (!is.na(pore_date)) {
        break
      }

      # Chercher une date incluse dans n'importe quelle cellule texte,
      # par exemple "Data_21_09-04-26_aHL_...".
      for (current_col in names(candidate_data)) {

        if (
          is.character(candidate_data[[current_col]]) ||
          is.factor(candidate_data[[current_col]])
        ) {

          candidate_text <- paste(
            as.character(candidate_data[[current_col]]),
            collapse = " "
          )

          candidate_date <- extract_date_from_text(
            candidate_text
          )

          if (!is.na(candidate_date)) {
            pore_date <- candidate_date
            break
          }
        }
      }

      if (!is.na(pore_date)) {
        break
      }
    }
  }

  # 4) Dernier recours : chemin / nom du fichier.
  if (is.na(pore_date)) {
    pore_date <- extract_date_from_text(
      file_path
    )
  }

  tibble::tibble(
    pore_id = pore_name,
    pore_duration_s = duration_s,
    Pore_duration_min = ifelse(
      is.finite(duration_s),
      duration_s / 60,
      NA_real_
    ),

    # Ceci est EXACTEMENT la baseline horizontale
    # utilisee sur le panel trace de la figure.
    Baseline_pA = abs(
      baseline_for_plot
    ),

    Baseline_min_pA = ifelse(
      is.finite(baseline_min),
      abs(baseline_min),
      NA_real_
    ),

    Baseline_max_pA = ifelse(
      is.finite(baseline_max),
      abs(baseline_max),
      NA_real_
    ),

    Lipid_concentration_mg_mL =
      lipid_concentration,

    Date = pore_date
  )
}


# ============================================================
# 7. BUILD PORE IDS EXACTLY LIKE FIGURE
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
        sub(
          "^P",
          "",
          pore_id
        )
      )
    )
  ) |>
  dplyr::arrange(
    pore_number,
    pore_id
  ) |>
  # EXACTEMENT comme le script de figure :
  # si deux fichiers donnent le meme P#, seul le premier est garde.
  dplyr::distinct(
    pore_id,
    .keep_all = TRUE
  )

pore_names <- pore_info$pore_id


# ============================================================
# 8. EXACT EVENT FILTER USED BY FIGURE
# ============================================================

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
  )

event_counts <- tibble::tibble(
  pore_id = pore_names
) |>
  dplyr::left_join(
    events_plot |>
      dplyr::group_by(
        pore_id
      ) |>
      dplyr::summarise(
        Number_of_events =
          dplyr::n(),
        .groups = "drop"
      ),
    by = "pore_id"
  ) |>
  dplyr::mutate(
    Number_of_events =
      dplyr::coalesce(
        Number_of_events,
        0L
      )
  )

# EXACTEMENT la meme regle que la figure : >= 100.
kept_pores <- event_counts |>
  dplyr::filter(
    Number_of_events >=
      MIN_EVENTS_PER_PORE
  ) |>
  dplyr::pull(
    pore_id
  )

if (length(kept_pores) == 0) {
  stop(
    "Aucun pore ne passe le filtre de ",
    MIN_EVENTS_PER_PORE,
    " evenements."
  )
}

pore_info_kept <- pore_info |>
  dplyr::filter(
    pore_id %in% kept_pores
  )

event_counts_kept <- event_counts |>
  dplyr::filter(
    pore_id %in% kept_pores
  )


# ============================================================
# 9. TRACE / BASELINE / DURATION FOR THE SAME KEPT PORES
# ============================================================

trace_summary <- purrr::map2_dfr(
  pore_info_kept$file_path,
  pore_info_kept$pore_id,
  read_trace_summary
)


# ============================================================
# 10. FINAL TABLE
# ============================================================

final_table <- pore_info_kept |>
  dplyr::select(
    pore_id,
    pore_number,
    file_name
  ) |>
  dplyr::left_join(
    event_counts_kept,
    by = "pore_id"
  ) |>
  dplyr::left_join(
    trace_summary,
    by = "pore_id"
  ) |>
  dplyr::arrange(
    pore_number,
    pore_id
  ) |>
  dplyr::mutate(
    # Row 1, Row 2, Row 3... = vertical order in the final figure.
    Figure_row_number =
      dplyr::row_number(),

    Event_frequency_Hz =
      dplyr::if_else(
        is.finite(pore_duration_s) &
          pore_duration_s > 0,
        Number_of_events /
          pore_duration_s,
        NA_real_
      ),

    Baseline_range_pA =
      dplyr::case_when(
        is.finite(Baseline_min_pA) &
          is.finite(Baseline_max_pA) ~
          paste0(
            format(
              round(
                Baseline_min_pA,
                3
              ),
              nsmall = 3,
              trim = TRUE
            ),
            " - ",
            format(
              round(
                Baseline_max_pA,
                3
              ),
              nsmall = 3,
              trim = TRUE
            ),
            " pA"
          ),
        TRUE ~ NA_character_
      )
  ) |>
  dplyr::transmute(
    Peptide = peptide_id,

    # EXACT label written in the scatter title: P1, P2, P3...
    Pore_on_figure = pore_id,

    # Position verticale du pore dans la figure.
    Figure_row_number =
      Figure_row_number,

    Number_of_events =
      Number_of_events,

    Event_frequency_Hz =
      round(
        Event_frequency_Hz,
        5
      ),

    Pore_duration_min =
      round(
        Pore_duration_min,
        3
      ),

    # EXACT dashed baseline used by the figure's raw-current panel.
    Baseline_pA =
      round(
        Baseline_pA,
        3
      ),

    # Range WITHIN THIS PORE.
    Baseline_min_pA =
      round(
        Baseline_min_pA,
        3
      ),

    Baseline_max_pA =
      round(
        Baseline_max_pA,
        3
      ),

    Baseline_range_pA =
      Baseline_range_pA,

    Lipid_concentration_mg_mL =
      Lipid_concentration_mg_mL,

    Date =
      Date,

    Source_file =
      file_name
  )


# ============================================================
# 11. EXPORT EXCEL
# ============================================================

output_file <- file.path(
  output_folder,
  paste0(
    peptide_id,
    "_figure_matching_table_WITH_DATE.xlsx"
  )
)

wb_out <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  wb_out,
  "Figure_matching_table"
)

openxlsx::writeDataTable(
  wb_out,
  "Figure_matching_table",
  final_table,
  tableStyle = "TableStyleMedium2"
)

openxlsx::freezePane(
  wb_out,
  "Figure_matching_table",
  firstRow = TRUE
)

openxlsx::setColWidths(
  wb_out,
  "Figure_matching_table",
  cols = seq_len(
    ncol(final_table)
  ),
  widths = 18
)

openxlsx::setColWidths(
  wb_out,
  "Figure_matching_table",
  cols = which(
    names(final_table) ==
      "Source_file"
  ),
  widths = 40
)

if (
  "Date" %in%
    names(final_table) &&
  nrow(final_table) > 0
) {

  openxlsx::addStyle(
    wb_out,
    "Figure_matching_table",
    style = openxlsx::createStyle(
      numFmt = "dd/mm/yyyy"
    ),
    rows = 2:(
      nrow(final_table) + 1
    ),
    cols = which(
      names(final_table) ==
        "Date"
    ),
    gridExpand = TRUE,
    stack = TRUE
  )
}

openxlsx::saveWorkbook(
  wb_out,
  output_file,
  overwrite = TRUE
)


# ============================================================
# 12. CONSOLE CHECK
# ============================================================

message("\n==============================================")
message("TABLE CORRESPONDANT EXACTEMENT A LA FIGURE")
message("Peptide : ", peptide_id)
message(
  "Pores de la figure : ",
  paste(
    final_table$Pore_on_figure,
    collapse = ", "
  )
)
message(
  "Nombre de pores : ",
  nrow(final_table)
)
message(
  "Filtre exact : dwell 0.2-100 ms ; I/I0 0-100 % ; n >= 100"
)
message(
  "Fichier cree : ",
  output_file
)
message("==============================================\n")

print(
  final_table |>
    dplyr::select(
      Pore_on_figure,
      Figure_row_number,
      Number_of_events,
      Event_frequency_Hz,
      Baseline_pA,
      Baseline_range_pA
    )
)
