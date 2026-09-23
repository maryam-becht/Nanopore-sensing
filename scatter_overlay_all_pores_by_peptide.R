
DATA_DIR <- "C:/Users/mbech/Documents/Peptides"


OUTPUT_DIR <- file.path(DATA_DIR, "Scatter_overlay_all_pores")


MIN_DWELL_MS <- 0.5
MAX_DWELL_MS <- 100
MIN_I_OVER_I0 <- 0.00
MAX_I_OVER_I0 <- 1.00  # mettre 0.70 pour une figure limitee a 70 %


POINT_SIZE <- 0.55
POINT_ALPHA <- 0.40
SHOW_PEPTIDE_TITLE <- FALSE
PNG_DPI <- 300

MAX_POINTS_PER_PORE <- Inf


DEFAULT_DWELL_UNIT <- "ms"  # choix possibles : "ms", "s", "us"



required_packages <- c("readxl", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Packages manquants : ", paste(missing_packages, collapse = ", "),
      "\nInstalle-les une seule fois avec :\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
    ),
    call. = FALSE
  )
}



clean_colnames <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x[is.na(x)] <- ""
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  make.unique(x, sep = "_")
}

as_numeric_safe <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  
  x <- trimws(as.character(x))
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("%", "", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

first_matching_col <- function(names_vector, exact = character(), regex = character()) {
  exact_hit <- exact[exact %in% names_vector]
  if (length(exact_hit) > 0) return(exact_hit[1])
  
  for (pattern in regex) {
    regex_hit <- grep(pattern, names_vector, value = TRUE, perl = TRUE)
    if (length(regex_hit) > 0) return(regex_hit[1])
  }
  
  NA_character_
}

detect_columns <- function(column_names) {
  dwell_col <- first_matching_col(
    column_names,
    exact = c(
      "dwell_ms", "dwell_time_ms", "duration_ms", "event_duration_ms",
      "dwell", "dwell_time", "duration", "event_duration",
      "dwell_s", "dwell_time_s", "duration_s",
      "dwell_us", "dwell_time_us", "duration_us"
    ),
    regex = c("dwell.*(ms|millisecond)", "duration.*(ms|millisecond)",
              "dwell", "duration")
  )
  
  ratio_col <- first_matching_col(
    column_names,
    exact = c(
      "i_over_i0", "i_over_i_0", "i_i0", "i_i_0", "ioveri0",
      "i0_ratio", "current_ratio", "residual_current",
      "relative_current", "fractional_current", "normalized_current",
      "i_over_i0_percent", "i_over_i0_pct"
    ),
    regex = c(
      "^i.*over.*i_?0", "^i.*i_?0", "ioveri0",
      "residual.*current", "relative.*current", "normalized.*current"
    )
  )
  
  pore_col <- first_matching_col(
    column_names,
    exact = c(
      "pore_id", "pore", "pore_number", "pore_no", "pore_index",
      "physical_pore", "true_pore", "run_id", "recording_id", "file_id"
    ),
    regex = c("^pore(_|$)", "physical.*pore", "true.*pore", "^run_id$")
  )
  
  peptide_col <- first_matching_col(
    column_names,
    exact = c("peptide", "peptide_name", "analyte", "sample"),
    regex = c("peptide", "analyte")
  )
  
  event_current_col <- first_matching_col(
    column_names,
    exact = c(
      "event_current", "event_current_pa", "mean_event_current",
      "mean_current", "avg_current", "average_current", "amplitude"
    ),
    regex = c("event.*current", "avg.*current", "mean.*current")
  )
  
  baseline_col <- first_matching_col(
    column_names,
    exact = c(
      "i0", "i_0", "open_pore_current", "open_pore_current_pa",
      "baseline", "baseline_pa", "baseline_current", "opc"
    ),
    regex = c("open.*pore.*current", "baseline.*current", "^i_?0($|_)")
  )
  
  depth_col <- first_matching_col(
    column_names,
    exact = c("depth_pa", "blockade_depth_pa", "current_drop_pa", "delta_i"),
    regex = c("depth.*pa", "blockade.*depth", "current.*drop", "delta.*i")
  )
  
  list(
    dwell = dwell_col,
    ratio = ratio_col,
    pore = pore_col,
    peptide = peptide_col,
    event_current = event_current_col,
    baseline = baseline_col,
    depth = depth_col
  )
}

convert_dwell_to_ms <- function(values, column_name) {
  values <- as_numeric_safe(values)
  
  if (grepl("(^|_)(us|microsecond|microseconds)($|_)", column_name)) {
    return(values / 1000)
  }
  
  if (grepl("(^|_)(ms|millisecond|milliseconds)($|_)", column_name)) {
    return(values)
  }
  
  if (column_name %in% c("duration", "duration_s", "dwell_s", "dwell_time_s",
                         "event_duration_s") ||
      grepl("(^|_)(s|sec|second|seconds)$", column_name)) {
    return(values * 1000)
  }
  
  switch(
    DEFAULT_DWELL_UNIT,
    "s" = values * 1000,
    "us" = values / 1000,
    values
  )
}

infer_peptide_from_path <- function(path) {
  path_abs <- normalizePath(path, winslash = "/", mustWork = FALSE)
  data_abs <- normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE)
  parent <- basename(dirname(path_abs))
  stem <- tools::file_path_sans_ext(basename(path_abs))
  
  generic_folders <- c(
    tolower(basename(data_abs)), "data", "excel", "excels", "events",
    "results", "output", "outputs", "recording", "recordings", "r_data"
  )
  
  label <- if (tolower(parent) %in% generic_folders || dirname(path_abs) == data_abs) {
    stem
  } else {
    parent
  }
  
  label <- gsub(
    "(?i)(_nanopore)?_r_data.*$|_all_pores.*$|_event(s)?(_data)?$|_analysis$",
    "", label, perl = TRUE
  )
  # Permet aussi de regrouper asynY125_pore1.xlsx, asynY125_pore2.xlsx, etc.
  label <- gsub("(?i)[_-](pore|p)[_-]*[0-9]+.*$", "", label, perl = TRUE)
  label <- gsub("[_-]+$", "", label)
  
  if (!nzchar(label)) stem else label
}

infer_pore_from_source <- function(path, sheet_name, total_sheets) {
  # Chaque fichier Excel correspond maintenant a un pore physique.
  # Exemple : asynpS87_P3_nanopore_R_data.xlsx -> P3
  stem <- tools::file_path_sans_ext(basename(path))
  pore_match <- regexpr(
    "(?i)(^|[_-])p[0-9]+(?=([_-]|$))",
    stem,
    perl = TRUE
  )
  
  if (pore_match[1] != -1) {
    pore_label <- regmatches(stem, pore_match)
    pore_label <- gsub("^[_-]+", "", pore_label)
    return(toupper(pore_label))
  }
  
  stem
}

calculate_i_over_i0 <- function(data, columns) {
  if (!is.na(columns$ratio)) {
    ratio <- as_numeric_safe(data[[columns$ratio]])
    finite_ratio <- ratio[is.finite(ratio)]
    
    # Les fichiers peuvent stocker I/I0 soit entre 0 et 1, soit entre 0 et 100.
    if (length(finite_ratio) > 0 && stats::median(abs(finite_ratio), na.rm = TRUE) > 1.5) {
      ratio <- ratio / 100
    }
    return(ratio)
  }
  
  if (!is.na(columns$event_current) && !is.na(columns$baseline)) {
    event_current <- abs(as_numeric_safe(data[[columns$event_current]]))
    baseline <- abs(as_numeric_safe(data[[columns$baseline]]))
    baseline[baseline == 0] <- NA_real_
    return(event_current / baseline)
  }
  
  if (!is.na(columns$depth) && !is.na(columns$baseline)) {
    depth <- abs(as_numeric_safe(data[[columns$depth]]))
    baseline <- abs(as_numeric_safe(data[[columns$baseline]]))
    baseline[baseline == 0] <- NA_real_
    return(1 - depth / baseline)
  }
  
  rep(NA_real_, nrow(data))
}

sheet_is_usable <- function(columns) {
  has_dwell <- !is.na(columns$dwell)
  has_ratio <- !is.na(columns$ratio)
  can_compute_from_current <- !is.na(columns$event_current) && !is.na(columns$baseline)
  can_compute_from_depth <- !is.na(columns$depth) && !is.na(columns$baseline)
  
  has_dwell && (has_ratio || can_compute_from_current || can_compute_from_depth)
}

read_one_sheet <- function(path, sheet_name, total_sheets) {
  header_preview <- tryCatch(
    readxl::read_excel(path, sheet = sheet_name, n_max = 5),
    error = function(e) NULL
  )
  
  if (is.null(header_preview) || ncol(header_preview) == 0) return(NULL)
  
  names(header_preview) <- clean_colnames(names(header_preview))
  preview_columns <- detect_columns(names(header_preview))
  
  # Evite de charger integralement les sheets de trace ou de resume inutiles.
  if (!sheet_is_usable(preview_columns)) return(NULL)
  
  raw <- tryCatch(
    readxl::read_excel(path, sheet = sheet_name),
    error = function(e) {
      warning("Sheet ignoree : ", basename(path), " / ", sheet_name,
              " — ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  
  names(raw) <- clean_colnames(names(raw))
  columns <- detect_columns(names(raw))
  if (!sheet_is_usable(columns)) return(NULL)
  
  dwell_ms <- convert_dwell_to_ms(raw[[columns$dwell]], columns$dwell)
  i_over_i0 <- calculate_i_over_i0(raw, columns)
  
  inferred_peptide <- infer_peptide_from_path(path)
  peptide <- if (!is.na(columns$peptide)) {
    as.character(raw[[columns$peptide]])
  } else {
    rep(inferred_peptide, nrow(raw))
  }
  peptide[is.na(peptide) | !nzchar(trimws(peptide))] <- inferred_peptide
  
  # Important : un fichier Excel = un pore. On ignore donc pore_id/run_id
  # a l'interieur du fichier, qui peuvent correspondre a des segments ou runs.
  inferred_pore <- infer_pore_from_source(path, sheet_name, total_sheets)
  pore <- rep(inferred_pore, nrow(raw))
  
  data.frame(
    peptide = trimws(peptide),
    pore_raw = trimws(pore),
    dwell_ms = dwell_ms,
    i_over_i0 = i_over_i0,
    source_file = basename(path),
    source_sheet = sheet_name,
    stringsAsFactors = FALSE
  )
}

make_safe_filename <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x[is.na(x)] <- "peptide"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "peptide")
}

natural_order <- function(x) {
  has_number <- grepl("[0-9]+", x)
  first_number <- rep(Inf, length(x))
  first_number[has_number] <- suppressWarnings(as.numeric(
    sub("^[^0-9]*([0-9]+).*$", "\\1", x[has_number])
  ))
  first_number[is.na(first_number)] <- Inf
  order(first_number, tolower(x))
}

standardize_pore_labels <- function(data) {
  original_pores <- unique(data$pore_raw)
  original_pores <- original_pores[natural_order(original_pores)]
  
  already_simple <- grepl(
    "^(pore|p)?[ _-]*[0-9]+$",
    original_pores,
    ignore.case = TRUE
  )
  labels <- if (all(already_simple)) {
    numbers <- suppressWarnings(as.integer(gsub("[^0-9]", "", original_pores)))
    paste0("P", numbers)
  } else {
    paste0("P", seq_along(original_pores))
  }
  
  if (anyDuplicated(labels)) labels <- paste0("P", seq_along(original_pores))
  
  mapping <- stats::setNames(labels, original_pores)
  data$pore <- unname(mapping[data$pore_raw])
  data$pore <- factor(data$pore, levels = labels)
  data
}

sample_for_display <- function(data) {
  if (is.infinite(MAX_POINTS_PER_PORE)) return(data)
  
  set.seed(125)
  groups <- split(data, interaction(data$peptide, data$pore, drop = TRUE))
  sampled <- lapply(groups, function(group) {
    if (nrow(group) <= MAX_POINTS_PER_PORE) return(group)
    group[sample.int(nrow(group), MAX_POINTS_PER_PORE), , drop = FALSE]
  })
  
  do.call(rbind, sampled)
}

make_overlay_plot <- function(data, peptide_name) {
  data <- standardize_pore_labels(data)
  plot_data <- sample_for_display(data)
  
  pore_levels <- levels(data$pore)

  article_palette <- c(
    "#4E79A7",  # bleu doux
    "#E07B62",  # corail doux
    "#59A14F",  # vert assourdi
    "#9C6FAE",  # violet doux
    "#6B6B6B",  # gris neutre
    "#76A5AF",  # bleu-gris
    "#B07A5A",  # terracotta
    "#7E6E8E"   # gris-violet
  )
  colours <- rep(article_palette, length.out = length(pore_levels))
  names(colours) <- pore_levels
  
  dwell_breaks <- c(0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
  dwell_breaks <- dwell_breaks[
    dwell_breaks >= MIN_DWELL_MS & dwell_breaks <= MAX_DWELL_MS
  ]
  
  iio_step <- if (MAX_I_OVER_I0 <= 0.70) 0.10 else 0.20
  iio_breaks <- seq(MIN_I_OVER_I0, MAX_I_OVER_I0, by = iio_step)
  
  plot_title <- if (SHOW_PEPTIDE_TITLE) peptide_name else NULL
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = dwell_ms, y = i_over_i0, colour = pore)
  ) +
    ggplot2::geom_point(
      size = POINT_SIZE,
      alpha = POINT_ALPHA,
      shape = 16,
      stroke = 0,
      na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(values = colours, name = "Pore", drop = FALSE) +
    ggplot2::scale_x_log10(
      limits = c(MIN_DWELL_MS, MAX_DWELL_MS),
      breaks = dwell_breaks,
      labels = format(dwell_breaks, trim = TRUE, scientific = FALSE),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(MIN_I_OVER_I0, MAX_I_OVER_I0),
      breaks = iio_breaks,
      labels = function(x) format(round(100 * x), trim = TRUE, scientific = FALSE),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Dwell time (ms)",
      y = expression(I/I[0] ~ "(%)")
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = grid::unit(1.8, "mm"),
      axis.text = ggplot2::element_text(colour = "black", size = 9),
      axis.title = ggplot2::element_text(colour = "black", size = 11),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      legend.key = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}


# ---------------------------
# 4. LECTURE DES EXCEL
# ---------------------------

if (!dir.exists(DATA_DIR)) {
  stop(
    paste0(
      "DATA_DIR n'existe pas :\n", DATA_DIR,
      "\n\nModifie DATA_DIR en haut du script avec le chemin de ton dossier Excel."
    ),
    call. = FALSE
  )
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

peptide_dirs <- list.dirs(
  DATA_DIR,
  recursive = FALSE,
  full.names = TRUE
)

output_abs <- normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE)
peptide_dirs <- peptide_dirs[
  normalizePath(peptide_dirs, winslash = "/", mustWork = FALSE) != output_abs
]
peptide_dirs <- peptide_dirs[order(tolower(basename(peptide_dirs)))]

if (length(peptide_dirs) == 0) {
  stop(
    paste0(
      "Aucun sous-dossier de peptide trouve dans :\n", DATA_DIR,
      "\nLa structure attendue est Peptides/nom_du_peptide/fichiers_P1_P2_etc."
    ),
    call. = FALSE
  )
}

summary_parts <- list()
n_figures <- 0L

for (peptide_dir in peptide_dirs) {
  peptide_name <- basename(peptide_dir)
  
  excel_files <- list.files(
    peptide_dir,
    pattern = "\\.(xlsx|xls)$",
    recursive = FALSE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  excel_files <- excel_files[!grepl("(^|[/\\\\])~\\$", excel_files)]
  excel_files <- excel_files[order(tolower(basename(excel_files)))]
  
  if (length(excel_files) == 0) next
  
  cat(
    "\nPeptide :", peptide_name,
    "|", length(excel_files), "fichier(s) Excel =",
    length(excel_files), "pore(s)\n"
  )
  
  peptide_parts <- list()
  part_index <- 0L
  
  for (file_path in excel_files) {
    sheet_names <- tryCatch(
      readxl::excel_sheets(file_path),
      error = function(e) {
        warning("Fichier ignore : ", basename(file_path), " — ",
                conditionMessage(e), call. = FALSE)
        character()
      }
    )
    
    if (length(sheet_names) == 0) next
    
    for (sheet_name in sheet_names) {
      part <- read_one_sheet(file_path, sheet_name, length(sheet_names))
      if (!is.null(part) && nrow(part) > 0) {
        part$peptide <- peptide_name
        part_index <- part_index + 1L
        peptide_parts[[part_index]] <- part
        cat(
          "  Lu :", basename(file_path), "|", sheet_name,
          "|", format(nrow(part), big.mark = ","), "evenements\n"
        )
      }
    }
  }
  
  if (length(peptide_parts) == 0) {
    warning(
      "Aucune sheet d'evenements reconnue pour ", peptide_name,
      ". Colonnes attendues : dwell_ms + I_over_I0, ou duration + avg_current + opc.",
      call. = FALSE
    )
    next
  }
  
  peptide_raw <- do.call(rbind, peptide_parts)
  rm(peptide_parts)
  invisible(gc(verbose = FALSE))
  peptide_raw <- standardize_pore_labels(peptide_raw)
  
  peptide_data <- peptide_raw[
    is.finite(peptide_raw$dwell_ms) &
      is.finite(peptide_raw$i_over_i0) &
      peptide_raw$dwell_ms >= MIN_DWELL_MS &
      peptide_raw$dwell_ms <= MAX_DWELL_MS &
      peptide_raw$i_over_i0 >= MIN_I_OVER_I0 &
      peptide_raw$i_over_i0 <= MAX_I_OVER_I0,
    ,
    drop = FALSE
  ]
  
  if (nrow(peptide_data) == 0) {
    warning(
      "Aucun evenement entre ", MIN_DWELL_MS, " et ", MAX_DWELL_MS,
      " ms pour ", peptide_name, ".",
      call. = FALSE
    )
    next
  }
  
  counts_in_excel <- aggregate(
    dwell_ms ~ peptide + pore,
    data = peptide_raw,
    FUN = length
  )
  names(counts_in_excel)[names(counts_in_excel) == "dwell_ms"] <-
    "n_events_in_excel"
  
  counts_plotted <- aggregate(
    dwell_ms ~ peptide + pore,
    data = peptide_data,
    FUN = length
  )
  names(counts_plotted)[names(counts_plotted) == "dwell_ms"] <-
    "n_events_plotted"
  
  event_counts <- merge(
    counts_in_excel,
    counts_plotted,
    by = c("peptide", "pore"),
    all.x = TRUE,
    sort = FALSE
  )
  event_counts$n_events_plotted[is.na(event_counts$n_events_plotted)] <- 0L
  summary_parts[[length(summary_parts) + 1L]] <- event_counts
  
  rm(peptide_raw)
  invisible(gc(verbose = FALSE))
  
  figure <- make_overlay_plot(peptide_data, peptide_name)
  safe_name <- make_safe_filename(peptide_name)
  
  png_path <- file.path(
    OUTPUT_DIR,
    paste0(safe_name, "_dwell_IoverI0_overlay_all_pores.png")
  )
  
  ggplot2::ggsave(
    filename = png_path,
    plot = figure,
    width = 6.6,
    height = 5.2,
    units = "in",
    dpi = PNG_DPI,
    bg = "white"
  )
  
  n_figures <- n_figures + 1L
  cat(
    "Figure creee :", peptide_name,
    "|", length(excel_files), "pores dans le dossier",
    "|", length(unique(peptide_data$pore)), "pores affiches",
    "|", format(nrow(peptide_data), big.mark = ","), "evenements\n"
  )
  
  rm(peptide_data, figure)
  if (exists("part", inherits = FALSE)) rm(part)
  invisible(gc(verbose = FALSE))
}

if (n_figures == 0) {
  stop(
    "Aucune figure n'a ete creee. Verifie les colonnes des fichiers Excel.",
    call. = FALSE
  )
}

event_count_summary <- do.call(rbind, summary_parts)
utils::write.csv(
  event_count_summary,
  file.path(OUTPUT_DIR, "event_counts_by_peptide_and_pore.csv"),
  row.names = FALSE
)

cat(
  "\nTermine :", n_figures, "figure(s) PNG creee(s).\n",
  "Resultats enregistres dans :\n", OUTPUT_DIR, "\n",
  sep = ""
)
