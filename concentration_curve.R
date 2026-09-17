# ============================================================
# NANOpore CONCENTRATION CURVE
# AT8 P1
# 200 pM -> 1 nM -> 50 nM -> 100 nM -> 500 nM -> 1 uM
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)


# ============================================================
# 1. FILE
# ============================================================

FILE <- rFILE <- r"(C:\Users\mbech\Documents\Peptides\taupThr217\taupThr217_P1_nanopore_R_data.xlsx)"

OUTPUT_DIR <- dirname(FILE)


# ============================================================
# 2. READ EVENTS
# ============================================================

events <- read_excel(
  FILE,
  sheet = "events"
)


events <- events %>%
  mutate(
    event_time_s = as.numeric(start_s_global)
  ) %>%
  filter(
    !is.na(event_time_s)
  )


# ============================================================
# 3. CHECK TIME RANGE
# ============================================================

cat("\nFirst event:",
    min(events$event_time_s, na.rm = TRUE),
    "s\n")

cat("Last event:",
    max(events$event_time_s, na.rm = TRUE),
    "s\n")


# ============================================================
# 4. CONCENTRATION WINDOWS
# ============================================================
#
# APPROXIMATION:
# recording ~300 s
# 6 concentrations
# ~50 s per concentration
#
# CHANGE THESE TIMES if you remember the real additions.
# ============================================================

windows <- data.frame(
  
  concentration = c(
    "200 pM",
    "1 nM",
    "50 nM",
    "100 nM",
    "500 nM",
    "1 µM"
  ),
  
  concentration_nM = c(
    0.2,
    1,
    50,
    100,
    500,
    1000
  ),
  
  start_s = c(
    0,
    50,
    100,
    150,
    200,
    250
  ),
  
  end_s = c(
    50,
    100,
    150,
    200,
    250,
    300
  )
)


# ============================================================
# 5. COUNT EVENTS
# ============================================================

results <- windows %>%
  rowwise() %>%
  mutate(
    
    N_events = sum(
      events$event_time_s >= start_s &
        events$event_time_s < end_s,
      na.rm = TRUE
    ),
    
    duration_s = end_s - start_s,
    
    frequency_Hz = N_events / duration_s
    
  ) %>%
  ungroup()


# ============================================================
# 6. PRINT RESULTS
# ============================================================

print(results)


# ============================================================
# 7. SAVE RESULTS
# ============================================================

write_xlsx(
  results,
  file.path(
    OUTPUT_DIR,
    "AT8_P1_concentration_curve_results.xlsx"
  )
)


# ============================================================
# 8. FIGURE
# ============================================================

plot_data <- results %>%
  filter(
    frequency_Hz > 0
  )


p <- ggplot(
  plot_data,
  aes(
    x = concentration_nM,
    y = frequency_Hz
  )
) +
  
  geom_line(
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey40"
  ) +
  
  geom_point(
    shape = 21,
    size = 4,
    fill = "white",
    colour = "black",
    stroke = 1
  ) +
  
  scale_x_log10(
    breaks = c(
      0.2,
      1,
      50,
      100,
      500,
      1000
    ),
    
    labels = c(
      "0.2",
      "1",
      "50",
      "100",
      "500",
      "1000"
    )
  ) +
  
  scale_y_log10() +
  
  labs(
    x = "Peptide concentration (nM)",
    y = expression(
      f[sig]~"(s"^-1*")"
    )
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text = element_text(
      colour = "black"
    ),
    
    axis.title = element_text(
      colour = "black"
    ),
    
    panel.grid = element_blank()
  )


print(p)


# ============================================================
# 9. SAVE FIGURE
# ============================================================

ggsave(
  file.path(
    OUTPUT_DIR,
    "taupThr217_P1_concentration_curve.png"
  ),
  plot = p,
  width = 5.5,
  height = 4.5,
  dpi = 600
)
