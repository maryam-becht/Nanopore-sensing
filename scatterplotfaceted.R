
required_packages <- c(
  "openxlsx",
  "dplyr",
  "purrr",
  "ggplot2",
  "stringr",
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

library(openxlsx)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(tibble)


# ============================================================
# 1. SÉLECTION DU DOSSIER CONTENANT TOUS LES PEPTIDES
# ============================================================

if (
  interactive() &&
  .Platform$OS.type == "windows"
) {
  
  peptides_root_folder <- choose.dir(
    default = "C:/Users/mbech/Documents/Peptides",
    caption = paste(
      "Choisis le dossier contenant",
      "tous les dossiers peptides"
    )
  )
  
} else {
  
  peptides_root_folder <- readline(
    paste0(
      "Entre le chemin complet du dossier ",
      "contenant tous les peptides : "
    )
  )
}

if (
  is.na(peptides_root_folder) ||
  peptides_root_folder == ""
) {
  
  stop(
    "Aucun dossier sélectionné."
  )
}

peptides_root_folder <- normalizePath(
  peptides_root_folder,
  winslash = "/",
  mustWork = TRUE
)

message(
  "Dossier parent : ",
  peptides_root_folder
)


# ============================================================
# 2. PARAMÈTRES
# ============================================================

# Dwell time affiché
MIN_DWELL_MS <- 0.5
MAX_DWELL_MS <- 100

# I/I0 minimum
MIN_I_OVER_I0 <- 0

# Les deux limites à générer
I_OVER_I0_LIMITS <- c(
  100,
  70
)

# Nombre minimum d'événements par pore après les filtres
#
# Garder 1 pour retrouver tous les pores visibles,
# y compris ceux qui ont très peu d'événements.
MIN_EVENTS_PER_PORE <- 1

# Nombre maximal de colonnes dans les facettes
MAX_FACET_COLUMNS <- 4

# Points des figures facettées
FACET_POINT_COLOUR <- "grey10"
FACET_POINT_ALPHA <- 0.50
FACET_POINT_SIZE <- 0.75

# Points des overlays
OVERLAY_POINT_ALPHA <- 0.45
OVERLAY_POINT_SIZE <- 0.75

# Qualité d'export
PNG_DPI <- 600

# TRUE = PNG et PDF
# FALSE = PNG seulement
SAVE_PDF <- FALSE


# ============================================================
# 3. DOSSIER GÉNÉRAL DE SORTIE
# ============================================================

general_output_folder <- file.path(
  peptides_root_folder,
  "R_all_peptides_output"
)

dir.create(
  general_output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

message(
  "Dossier de sortie : ",
  general_output_folder
)


# ============================================================
# 4. FONCTIONS UTILITAIRES
# ============================================================


# ------------------------------------------------------------
# Nom de fichier compatible Windows
# ------------------------------------------------------------

make_safe_filename <- function(x) {
  
  x |>
    stringr::str_replace_all(
      "[^A-Za-z0-9_-]+",
      "_"
    ) |>
    stringr::str_replace_all(
      "_+",
      "_"
    ) |>
    stringr::str_remove(
      "^_"
    ) |>
    stringr::str_remove(
      "_$"
    )
}


# ------------------------------------------------------------
# Vérifier si un dossier contient directement des Excel
# ------------------------------------------------------------

folder_contains_excel <- function(folder_path) {
  
  excel_files <- list.files(
    path = folder_path,
    pattern = "\\.(xlsx|xlsm)$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  
  excel_files <- excel_files[
    !stringr::str_detect(
      basename(excel_files),
      "^~\\$"
    )
  ]
  
  length(excel_files) > 0
}


# ------------------------------------------------------------
# Identifier correctement le pore dans le nom du fichier
#
# Important :
# Cette fonction ne prend pas pS87, pY125 ou pThr217
# pour des numéros de pores.
#
# Exemple :
# asynpS87_P1_nanopore_R_data.xlsx -> P1
# ------------------------------------------------------------

extract_pore_id_from_filename <- function(file_path) {
  
  filename_without_extension <- tools::file_path_sans_ext(
    basename(file_path)
  )
  
  
  # ----------------------------------------------------------
  # Cas principal :
  # _P1_nanopore_R_data
  # ----------------------------------------------------------
  
  main_match <- stringr::str_match(
    filename_without_extension,
    stringr::regex(
      "_(P[0-9]+)_nanopore_R_data$",
      ignore_case = TRUE
    )
  )
  
  if (!is.na(main_match[1, 2])) {
    
    return(
      toupper(
        main_match[1, 2]
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Cas général :
  # P1 doit être un élément distinct entouré de "_" ou "-"
  #
  # Exemples :
  # peptide_P1_events
  # peptide-P1-events
  # ----------------------------------------------------------
  
  fallback_match <- stringr::str_match(
    filename_without_extension,
    stringr::regex(
      "(?:^|[_-])(P[0-9]+)(?=[_-]|$)",
      ignore_case = TRUE
    )
  )
  
  if (!is.na(fallback_match[1, 2])) {
    
    return(
      toupper(
        fallback_match[1, 2]
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Aucun pore détecté
  # ----------------------------------------------------------
  
  NA_character_
}


# ------------------------------------------------------------
# Numéro numérique du pore
#
# P1  -> 1
# P10 -> 10
# ------------------------------------------------------------

extract_pore_number <- function(pore_id) {
  
  suppressWarnings(
    as.integer(
      stringr::str_extract(
        pore_id,
        "[0-9]+"
      )
    )
  )
}


# ------------------------------------------------------------
# Lire la première ligne d'une feuille Excel
# ------------------------------------------------------------

get_excel_header <- function(
    workbook,
    sheet_name
) {
  
  header_data <- openxlsx::read.xlsx(
    workbook,
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


# ------------------------------------------------------------
# Lire les événements d'un fichier = un pore
#
# Seules les colonnes nécessaires sont lues :
#   dwell_total_ms
#   I_over_I0_pct
# ------------------------------------------------------------

read_events_one_pore <- function(
    file_path,
    peptide_name,
    pore_name
) {
  
  workbook <- openxlsx::loadWorkbook(
    file_path
  )
  
  sheet_names <- names(
    workbook
  )
  
  
  # Trouver la feuille events sans tenir compte
  # des majuscules et minuscules
  
  events_sheet <- sheet_names[
    tolower(
      trimws(
        sheet_names
      )
    ) == "events"
  ]
  
  if (length(events_sheet) == 0) {
    
    stop(
      "La feuille 'events' est absente de ",
      basename(file_path),
      ". Feuilles disponibles : ",
      paste(
        sheet_names,
        collapse = ", "
      )
    )
  }
  
  events_sheet <- events_sheet[1]
  
  
  # Lire les en-têtes
  
  header <- get_excel_header(
    workbook,
    events_sheet
  )
  
  dwell_column <- match(
    "dwell_total_ms",
    header
  )
  
  iio_column <- match(
    "I_over_I0_pct",
    header
  )
  
  if (
    is.na(dwell_column) ||
    is.na(iio_column)
  ) {
    
    stop(
      "Colonnes dwell_total_ms et/ou I_over_I0_pct ",
      "manquantes dans ",
      basename(file_path)
    )
  }
  
  
  # Lire seulement les deux colonnes utiles
  
  selected_columns <- sort(
    unique(
      c(
        dwell_column,
        iio_column
      )
    )
  )
  
  events_data <- openxlsx::read.xlsx(
    workbook,
    sheet = events_sheet,
    cols = selected_columns,
    colNames = TRUE,
    check.names = FALSE
  )
  
  
  # Nettoyage et standardisation
  
  events_data |>
    dplyr::transmute(
      
      peptide_id = peptide_name,
      
      pore_id = pore_name,
      
      source_file = basename(
        file_path
      ),
      
      dwell_total_ms = suppressWarnings(
        as.numeric(
          .data[["dwell_total_ms"]]
        )
      ),
      
      I_over_I0_pct = suppressWarnings(
        as.numeric(
          .data[["I_over_I0_pct"]]
        )
      )
    ) |>
    dplyr::filter(
      is.finite(
        dwell_total_ms
      ),
      
      is.finite(
        I_over_I0_pct
      )
    )
}


# ============================================================
# 5. THÈME DES FIGURES
# ============================================================

scatter_theme <- ggplot2::theme_classic(
  base_size = 10
) +
  
  ggplot2::theme(
    
    axis.text = ggplot2::element_text(
      colour = "black",
      size = 8.5
    ),
    
    axis.text.x = ggplot2::element_text(
      colour = "black",
      size = 8
    ),
    
    axis.title = ggplot2::element_text(
      colour = "black",
      size = 10.5
    ),
    
    axis.line = ggplot2::element_line(
      colour = "black",
      linewidth = 0.35
    ),
    
    axis.ticks = ggplot2::element_line(
      colour = "black",
      linewidth = 0.30
    ),
    
    strip.background = ggplot2::element_blank(),
    
    strip.text = ggplot2::element_text(
      colour = "black",
      size = 9.5
    ),
    
    panel.border = ggplot2::element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.45
    ),
    
    panel.spacing = grid::unit(
      0.8,
      "lines"
    ),
    
    plot.margin = ggplot2::margin(
      t = 8,
      r = 10,
      b = 8,
      l = 10
    )
  )


# ============================================================
# 6. SAUVEGARDE PNG ET PDF
# ============================================================

save_plot_files <- function(
    plot_object,
    output_prefix,
    output_folder,
    width,
    height
) {
  
  png_path <- file.path(
    output_folder,
    paste0(
      output_prefix,
      ".png"
    )
  )
  
  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = PNG_DPI,
    bg = "white",
    limitsize = FALSE
  )
  
  if (SAVE_PDF) {
    
    pdf_path <- file.path(
      output_folder,
      paste0(
        output_prefix,
        ".pdf"
      )
    )
    
    ggplot2::ggsave(
      filename = pdf_path,
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      limitsize = FALSE
    )
  }
}


# ============================================================
# 7. TRAITER UN DOSSIER PEPTIDE
# ============================================================

process_one_peptide_folder <- function(
    peptide_folder
) {
  
  peptide_folder <- normalizePath(
    peptide_folder,
    winslash = "/",
    mustWork = TRUE
  )
  
  peptide_id <- basename(
    peptide_folder
  )
  
  safe_peptide_id <- make_safe_filename(
    peptide_id
  )
  
  message("")
  message(
    "================================================"
  )
  message(
    "PEPTIDE : ",
    peptide_id
  )
  message(
    "================================================"
  )
  
  
  # ----------------------------------------------------------
  # Dossier de sortie du peptide
  # ----------------------------------------------------------
  
  peptide_output_folder <- file.path(
    general_output_folder,
    safe_peptide_id
  )
  
  dir.create(
    peptide_output_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Trouver les fichiers Excel
  # ----------------------------------------------------------
  
  excel_files <- list.files(
    path = peptide_folder,
    pattern = "\\.(xlsx|xlsm)$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  
  excel_files <- excel_files[
    !stringr::str_detect(
      basename(excel_files),
      "^~\\$"
    )
  ]
  
  excel_files <- sort(
    excel_files
  )
  
  if (length(excel_files) == 0) {
    
    warning(
      "Aucun fichier Excel trouvé pour ",
      peptide_id
    )
    
    return(
      NULL
    )
  }
  
  message(
    length(excel_files),
    " fichier(s) Excel trouvé(s)."
  )
  
  
  # ==========================================================
  # IDENTIFIER CORRECTEMENT LES PORES
  # ==========================================================
  
  detected_pore_ids <- purrr::map_chr(
    excel_files,
    extract_pore_id_from_filename
  )
  
  pore_detection_check <- tibble::tibble(
    file_name = basename(
      excel_files
    ),
    
    detected_pore_id = detected_pore_ids
  )
  
  message(
    "Correspondance fichiers / pores :"
  )
  
  print(
    pore_detection_check
  )
  
  
  # ----------------------------------------------------------
  # Vérifier les fichiers sans identifiant de pore
  # ----------------------------------------------------------
  
  if (any(is.na(detected_pore_ids))) {
    
    files_without_pore_id <- basename(
      excel_files[
        is.na(detected_pore_ids)
      ]
    )
    
    stop(
      paste0(
        "Impossible d'identifier le pore dans :\n",
        paste(
          files_without_pore_id,
          collapse = "\n"
        ),
        "\n\nNom attendu, par exemple : ",
        "asynpS87_P1_nanopore_R_data.xlsx"
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Vérifier les pores en double
  # ----------------------------------------------------------
  
  duplicate_pores <- unique(
    detected_pore_ids[
      duplicated(
        detected_pore_ids
      )
    ]
  )
  
  if (length(duplicate_pores) > 0) {
    
    duplicate_details <- pore_detection_check |>
      dplyr::filter(
        detected_pore_id %in% duplicate_pores
      )
    
    print(
      duplicate_details
    )
    
    stop(
      "Plusieurs fichiers ont le même identifiant de pore : ",
      paste(
        duplicate_pores,
        collapse = ", "
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Tableau des pores trié numériquement
  # ----------------------------------------------------------
  
  pore_info <- tibble::tibble(
    file_path = excel_files,
    file_name = basename(
      excel_files
    ),
    pore_id = detected_pore_ids
  ) |>
    dplyr::mutate(
      pore_number = purrr::map_int(
        pore_id,
        extract_pore_number
      )
    ) |>
    dplyr::arrange(
      pore_number,
      pore_id
    )
  
  
  # ----------------------------------------------------------
  # Lire chaque Excel
  #
  # Une erreur sur un fichier n'arrête pas tous les peptides
  # ----------------------------------------------------------
  
  read_attempts <- purrr::map2(
    pore_info$file_path,
    pore_info$pore_id,
    function(
    current_file,
    current_pore
    ) {
      
      tryCatch(
        {
          list(
            data = read_events_one_pore(
              file_path = current_file,
              peptide_name = peptide_id,
              pore_name = current_pore
            ),
            
            error = NULL
          )
        },
        
        error = function(e) {
          
          list(
            data = NULL,
            error = conditionMessage(
              e
            )
          )
        }
      )
    }
  )
  
  
  # ----------------------------------------------------------
  # Statut de lecture des fichiers
  # ----------------------------------------------------------
  
  pore_info <- pore_info |>
    dplyr::mutate(
      
      read_status = purrr::map_chr(
        read_attempts,
        function(x) {
          
          if (is.null(x$error)) {
            "OK"
          } else {
            "ERROR"
          }
        }
      ),
      
      read_error = purrr::map_chr(
        read_attempts,
        function(x) {
          
          if (is.null(x$error)) {
            ""
          } else {
            x$error
          }
        }
      )
    )
  
  write.csv(
    pore_info |>
      dplyr::select(
        pore_id,
        file_name,
        read_status,
        read_error
      ),
    
    file = file.path(
      peptide_output_folder,
      paste0(
        safe_peptide_id,
        "_pore_mapping.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Garder les fichiers lus correctement
  # ----------------------------------------------------------
  
  valid_indices <- which(
    pore_info$read_status == "OK"
  )
  
  if (length(valid_indices) == 0) {
    
    warning(
      "Aucun fichier lisible pour ",
      peptide_id
    )
    
    return(
      NULL
    )
  }
  
  events_raw <- dplyr::bind_rows(
    purrr::map(
      read_attempts[
        valid_indices
      ],
      "data"
    )
  )
  
  valid_pore_info <- pore_info[
    valid_indices,
    ,
    drop = FALSE
  ]
  
  pore_names <- valid_pore_info$pore_id
  
  events_raw <- events_raw |>
    dplyr::mutate(
      pore_id = factor(
        pore_id,
        levels = pore_names
      )
    )
  
  if (nrow(events_raw) == 0) {
    
    warning(
      "Aucun événement valide pour ",
      peptide_id
    )
    
    return(
      NULL
    )
  }
  
  
  # ----------------------------------------------------------
  # Palette :
  # la même couleur est attribuée au même pore dans les
  # versions 100 % et 70 %
  # ----------------------------------------------------------
  
  pore_colours <- grDevices::hcl.colors(
    n = length(
      pore_names
    ),
    palette = "Dark 3"
  )
  
  names(
    pore_colours
  ) <- pore_names
  
  
  # ----------------------------------------------------------
  # Nombre d'événements présents dans chaque Excel
  # avant les filtres de visualisation
  # ----------------------------------------------------------
  
  events_in_excel_counts <- events_raw |>
    dplyr::mutate(
      pore_id = as.character(
        pore_id
      )
    ) |>
    dplyr::count(
      pore_id,
      name = "n_events_in_excel"
    )
  
  
  # ==========================================================
  # CRÉER LES FIGURES 100 % ET 70 %
  # ==========================================================
  
  all_limit_summaries <- purrr::map_dfr(
    I_OVER_I0_LIMITS,
    function(
    current_iio_limit
    ) {
      
      message(
        "Création des figures I/I0 0–",
        current_iio_limit,
        " %..."
      )
      
      
      # ------------------------------------------------------
      # Filtrage des événements
      # ------------------------------------------------------
      
      events_plot <- events_raw |>
        dplyr::filter(
          
          dwell_total_ms >= MIN_DWELL_MS,
          dwell_total_ms <= MAX_DWELL_MS,
          
          I_over_I0_pct >= MIN_I_OVER_I0,
          I_over_I0_pct <= current_iio_limit
        )
      
      
      # ------------------------------------------------------
      # Compter les événements par pore
      # ------------------------------------------------------
      
      event_counts <- tibble::tibble(
        pore_id = pore_names
      ) |>
        dplyr::left_join(
          events_in_excel_counts,
          by = "pore_id"
        ) |>
        dplyr::left_join(
          events_plot |>
            dplyr::mutate(
              pore_id = as.character(
                pore_id
              )
            ) |>
            dplyr::count(
              pore_id,
              name = "n_events_plotted"
            ),
          
          by = "pore_id"
        ) |>
        dplyr::mutate(
          
          n_events_in_excel = dplyr::coalesce(
            n_events_in_excel,
            0L
          ),
          
          n_events_plotted = dplyr::coalesce(
            n_events_plotted,
            0L
          ),
          
          peptide_id = peptide_id,
          
          I_over_I0_max_pct = current_iio_limit
        )
      
      
      # ------------------------------------------------------
      # Garder les pores ayant suffisamment d'événements
      # ------------------------------------------------------
      
      kept_pores <- event_counts |>
        dplyr::filter(
          n_events_plotted >= MIN_EVENTS_PER_PORE
        ) |>
        dplyr::pull(
          pore_id
        )
      
      if (length(kept_pores) == 0) {
        
        warning(
          "Aucun pore de ",
          peptide_id,
          " ne passe les filtres à ",
          current_iio_limit,
          " %."
        )
        
        return(
          event_counts
        )
      }
      
      events_plot <- events_plot |>
        dplyr::filter(
          as.character(
            pore_id
          ) %in% kept_pores
        ) |>
        dplyr::mutate(
          pore_id = factor(
            as.character(
              pore_id
            ),
            levels = kept_pores
          )
        )
      
      plotted_event_counts <- event_counts |>
        dplyr::filter(
          pore_id %in% kept_pores
        )
      
      
      # ------------------------------------------------------
      # Sauvegarder les comptes
      # ------------------------------------------------------
      
      write.csv(
        event_counts,
        
        file = file.path(
          peptide_output_folder,
          paste0(
            safe_peptide_id,
            "_event_counts_dwell0.5ms_IoverI0_0to",
            current_iio_limit,
            "pct.csv"
          )
        ),
        
        row.names = FALSE
      )
      
      
      # ------------------------------------------------------
      # Graduations axe X
      # ------------------------------------------------------
      
      if (current_iio_limit == 100) {
        
        x_breaks <- seq(
          0,
          100,
          by = 20
        )
        
      } else {
        
        x_breaks <- seq(
          0,
          current_iio_limit,
          by = 10
        )
      }
      
      
      # ======================================================
      # A. FACETED SCATTER PLOT
      # ======================================================
      
      number_of_pores <- length(
        kept_pores
      )
      
      facet_columns <- min(
        MAX_FACET_COLUMNS,
        number_of_pores
      )
      
      facet_rows <- ceiling(
        number_of_pores /
          facet_columns
      )
      
      facet_width <- max(
        4.2,
        2.55 * facet_columns + 0.9
      )
      
      facet_height <- max(
        3.8,
        2.35 * facet_rows + 0.9
      )
      
      
      # Ajouter n dans le titre de chaque facette
      
      facet_label_values <- setNames(
        paste0(
          plotted_event_counts$pore_id,
          "\nn = ",
          plotted_event_counts$n_events_plotted
        ),
        plotted_event_counts$pore_id
      )
      
      
      faceted_plot <- ggplot2::ggplot(
        events_plot,
        ggplot2::aes(
          x = I_over_I0_pct,
          y = dwell_total_ms
        )
      ) +
        
        ggplot2::geom_point(
          colour = FACET_POINT_COLOUR,
          alpha = FACET_POINT_ALPHA,
          size = FACET_POINT_SIZE
        ) +
        
        ggplot2::scale_x_continuous(
          limits = c(
            MIN_I_OVER_I0,
            current_iio_limit
          ),
          
          breaks = x_breaks,
          
          expand = ggplot2::expansion(
            mult = c(
              0,
              0.02
            )
          )
        ) +
        
        ggplot2::scale_y_log10(
          limits = c(
            MIN_DWELL_MS,
            MAX_DWELL_MS
          ),
          
          breaks = c(
            0.5,
            1,
            2,
            5,
            10,
            20,
            50,
            100
          ),
          
          labels = c(
            "0.5",
            "1",
            "2",
            "5",
            "10",
            "20",
            "50",
            "100"
          ),
          
          expand = ggplot2::expansion(
            mult = c(
              0,
              0.02
            )
          )
        ) +
        
        ggplot2::facet_wrap(
          ~ pore_id,
          ncol = facet_columns,
          scales = "fixed",
          drop = TRUE,
          
          labeller = ggplot2::as_labeller(
            facet_label_values
          )
        ) +
        
        ggplot2::labs(
          x = expression(
            I/I[0]~" (%)"
          ),
          
          y = "Dwell time (ms)"
        ) +
        
        scatter_theme +
        
        ggplot2::theme(
          legend.position = "none"
        )
      
      
      faceted_output_prefix <- paste0(
        safe_peptide_id,
        "_faceted_scatter_dwell0.5to100ms_IoverI0_0to",
        current_iio_limit,
        "pct"
      )
      
      save_plot_files(
        plot_object = faceted_plot,
        output_prefix = faceted_output_prefix,
        output_folder = peptide_output_folder,
        width = facet_width,
        height = facet_height
      )
      
      
      # ======================================================
      # B. OVERLAY SCATTER PLOT
      # ======================================================
      
      overlay_colours <- pore_colours[
        kept_pores
      ]
      
      overlay_labels <- setNames(
        paste0(
          plotted_event_counts$pore_id,
          " (n = ",
          plotted_event_counts$n_events_plotted,
          ")"
        ),
        plotted_event_counts$pore_id
      )
      
      overlay_plot <- ggplot2::ggplot(
        events_plot,
        ggplot2::aes(
          x = I_over_I0_pct,
          y = dwell_total_ms,
          colour = pore_id
        )
      ) +
        
        ggplot2::geom_point(
          alpha = OVERLAY_POINT_ALPHA,
          size = OVERLAY_POINT_SIZE
        ) +
        
        ggplot2::scale_colour_manual(
          values = overlay_colours,
          breaks = kept_pores,
          labels = overlay_labels[
            kept_pores
          ],
          drop = FALSE
        ) +
        
        ggplot2::scale_x_continuous(
          limits = c(
            MIN_I_OVER_I0,
            current_iio_limit
          ),
          
          breaks = x_breaks,
          
          expand = ggplot2::expansion(
            mult = c(
              0,
              0.02
            )
          )
        ) +
        
        ggplot2::scale_y_log10(
          limits = c(
            MIN_DWELL_MS,
            MAX_DWELL_MS
          ),
          
          breaks = c(
            0.5,
            1,
            2,
            5,
            10,
            20,
            50,
            100
          ),
          
          labels = c(
            "0.5",
            "1",
            "2",
            "5",
            "10",
            "20",
            "50",
            "100"
          ),
          
          expand = ggplot2::expansion(
            mult = c(
              0,
              0.02
            )
          )
        ) +
        
        ggplot2::labs(
          x = expression(
            I/I[0]~" (%)"
          ),
          
          y = "Dwell time (ms)",
          
          colour = "Pore"
        ) +
        
        scatter_theme +
        
        ggplot2::theme(
          
          legend.position = "right",
          
          legend.title = ggplot2::element_text(
            colour = "black",
            size = 9,
            face = "bold"
          ),
          
          legend.text = ggplot2::element_text(
            colour = "black",
            size = 8
          ),
          
          legend.key.height = grid::unit(
            0.45,
            "cm"
          )
        ) +
        
        ggplot2::guides(
          colour = ggplot2::guide_legend(
            override.aes = list(
              alpha = 1,
              size = 2
            )
          )
        )
      
      
      overlay_width <- if (
        number_of_pores <= 5
      ) {
        7
      } else {
        7.8
      }
      
      overlay_height <- 5
      
      overlay_output_prefix <- paste0(
        safe_peptide_id,
        "_overlay_scatter_dwell0.5to100ms_IoverI0_0to",
        current_iio_limit,
        "pct"
      )
      
      save_plot_files(
        plot_object = overlay_plot,
        output_prefix = overlay_output_prefix,
        output_folder = peptide_output_folder,
        width = overlay_width,
        height = overlay_height
      )
      
      message(
        "Figures terminées : ",
        peptide_id,
        " — I/I0 ≤ ",
        current_iio_limit,
        " %."
      )
      
      event_counts
    }
  )
  
  
  # ----------------------------------------------------------
  # Résumé du peptide
  # ----------------------------------------------------------
  
  write.csv(
    all_limit_summaries,
    
    file = file.path(
      peptide_output_folder,
      paste0(
        safe_peptide_id,
        "_all_plot_event_counts.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  message(
    "Peptide terminé : ",
    peptide_id
  )
  
  all_limit_summaries
}


# ============================================================
# 8. TROUVER TOUS LES DOSSIERS PEPTIDES
# ============================================================

candidate_folders <- list.dirs(
  path = peptides_root_folder,
  full.names = TRUE,
  recursive = FALSE
)

# Ne pas relire le dossier de sortie
candidate_folders <- candidate_folders[
  basename(candidate_folders) !=
    basename(general_output_folder)
]

peptide_folders <- candidate_folders[
  vapply(
    candidate_folders,
    folder_contains_excel,
    logical(1)
  )
]


# Permet aussi de sélectionner directement un dossier peptide
if (
  length(peptide_folders) == 0 &&
  folder_contains_excel(
    peptides_root_folder
  )
) {
  
  peptide_folders <- peptides_root_folder
}

if (length(peptide_folders) == 0) {
  
  stop(
    paste0(
      "Aucun dossier peptide contenant directement ",
      "des fichiers Excel n'a été trouvé dans : ",
      peptides_root_folder
    )
  )
}

peptide_folders <- sort(
  peptide_folders
)

message("")
message(
  "Nombre de peptides trouvés : ",
  length(peptide_folders)
)

print(
  basename(
    peptide_folders
  )
)


# ============================================================
# 9. TRAITER TOUS LES PEPTIDES
# ============================================================

all_peptides_summaries <- list()
processing_log <- list()

for (
  peptide_index in seq_along(
    peptide_folders
  )
) {
  
  current_folder <- peptide_folders[
    peptide_index
  ]
  
  current_peptide <- basename(
    current_folder
  )
  
  current_result <- tryCatch(
    {
      
      peptide_summary <- process_one_peptide_folder(
        current_folder
      )
      
      processing_log[[
        current_peptide
      ]] <- tibble::tibble(
        peptide_id = current_peptide,
        status = "OK",
        error = ""
      )
      
      peptide_summary
    },
    
    error = function(e) {
      
      error_message <- conditionMessage(
        e
      )
      
      warning(
        "Erreur pour ",
        current_peptide,
        " : ",
        error_message
      )
      
      processing_log[[
        current_peptide
      ]] <- tibble::tibble(
        peptide_id = current_peptide,
        status = "ERROR",
        error = error_message
      )
      
      NULL
    }
  )
  
  if (!is.null(current_result)) {
    
    all_peptides_summaries[[
      current_peptide
    ]] <- current_result
  }
}


# ============================================================
# 10. RÉCAPITULATIF GLOBAL
# ============================================================

global_processing_log <- dplyr::bind_rows(
  processing_log
)

write.csv(
  global_processing_log,
  
  file = file.path(
    general_output_folder,
    "all_peptides_processing_log.csv"
  ),
  
  row.names = FALSE
)

if (length(all_peptides_summaries) > 0) {
  
  global_event_summary <- dplyr::bind_rows(
    all_peptides_summaries
  )
  
  write.csv(
    global_event_summary,
    
    file = file.path(
      general_output_folder,
      "all_peptides_event_counts.csv"
    ),
    
    row.names = FALSE
  )
}


# ============================================================
# 11. FIN
# ============================================================

message("")
message(
  "================================================"
)
message(
  "ANALYSE TERMINÉE"
)
message(
  "================================================"
)

message(
  "Résultats enregistrés dans :"
)

message(
  general_output_folder
)
