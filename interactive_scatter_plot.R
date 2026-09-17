# ============================================================
# NANOPORE EVENT EXPLORER
#
# X = I/I0 (%)
# Y = Dwell time (ms, log scale)
#
# Click on a scatter point -> display exact event shape
# ============================================================

# ============================================================
# PACKAGES
# ============================================================

library(shiny)
library(plotly)
library(dplyr)
library(readr)


# ============================================================
# SELECT EVENT SHAPES CSV
# ============================================================

cat("\nChoose the *_event_shapes.csv file...\n")

csv_file <- file.choose()

cat("\nReading:\n")
cat(csv_file, "\n\n")


# ============================================================
# READ DATA
# ============================================================

traces <- read_csv(
  csv_file,
  show_col_types = FALSE
)


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "event_id",
  "file_name",
  "channel",
  "sweep",
  "time_relative_ms",
  "current_pA",
  "baseline_pA",
  "I_over_I0_trace_pct",
  "inside_event",
  "dwell_total_ms",
  "event_I_over_I0_pct"
)

missing_columns <- setdiff(
  required_columns,
  names(traces)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing columns:\n",
      paste(missing_columns, collapse = "\n")
    )
  )
}


# ============================================================
# EVENT SUMMARY
# ============================================================

events <- traces |>
  distinct(
    event_id,
    file_name,
    channel,
    sweep,
    dwell_total_ms,
    event_I_over_I0_pct,
    baseline_pA
  ) |>
  arrange(event_id)

cat(nrow(events), "events loaded.\n")


# ============================================================
# CLEAN LIMITS FOR FILTERS
# ============================================================

DWELL_MIN_DATA <- min(events$dwell_total_ms, na.rm = TRUE)
DWELL_MAX_DATA <- max(events$dwell_total_ms, na.rm = TRUE)

RATIO_MIN_DATA <- min(events$event_I_over_I0_pct, na.rm = TRUE)
RATIO_MAX_DATA <- max(events$event_I_over_I0_pct, na.rm = TRUE)

# clean slider display values
DWELL_MIN_UI <- floor(DWELL_MIN_DATA * 1000) / 1000
DWELL_MAX_UI <- ceiling(DWELL_MAX_DATA * 1000) / 1000

RATIO_MIN_UI <- floor(RATIO_MIN_DATA)
RATIO_MAX_UI <- ceiling(RATIO_MAX_DATA)


# ============================================================
# Y-AXIS TICKS FOR LOG DWELL AXIS
# ============================================================

candidate_ticks <- c(
  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10
)

dwell_ticks <- candidate_ticks[
  candidate_ticks >= DWELL_MIN_UI &
    candidate_ticks <= DWELL_MAX_UI
]

if (!any(dwell_ticks == 1)) {
  dwell_ticks <- sort(unique(c(dwell_ticks, 1)))
}

dwell_tick_text <- as.character(dwell_ticks)


# ============================================================
# USER INTERFACE
# ============================================================

ui <- fluidPage(
  
  titlePanel("Nanopore Event Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      width = 3,
      
      h4("Filters"),
      
      sliderInput(
        inputId = "ratio_filter",
        label = "I/I0 (%)",
        min = RATIO_MIN_UI,
        max = RATIO_MAX_UI,
        value = c(RATIO_MIN_UI, RATIO_MAX_UI),
        step = 1
      ),
      
      sliderInput(
        inputId = "dwell_filter",
        label = "Dwell time (ms)",
        min = DWELL_MIN_UI,
        max = DWELL_MAX_UI,
        value = c(DWELL_MIN_UI, DWELL_MAX_UI),
        step = 0.001
      ),
      
      hr(),
      
      h4("Event navigation"),
      
      numericInput(
        inputId = "event_number",
        label = "Event ID",
        value = min(events$event_id),
        min = min(events$event_id),
        max = max(events$event_id),
        step = 1
      ),
      
      actionButton(
        inputId = "go_event",
        label = "Go to event"
      ),
      
      br(),
      br(),
      
      fluidRow(
        column(
          width = 6,
          actionButton(
            inputId = "previous_event",
            label = "← Previous"
          )
        ),
        column(
          width = 6,
          actionButton(
            inputId = "next_event",
            label = "Next →"
          )
        )
      ),
      
      hr(),
      
      h4("Trace display"),
      
      radioButtons(
        inputId = "trace_mode",
        label = NULL,
        choices = c(
          "I/I0 (%)" = "ratio",
          "Current (pA)" = "current"
        ),
        selected = "ratio"
      ),
      
      hr(),
      
      strong(
        textOutput("event_count")
      )
      
    ),
    
    mainPanel(
      
      width = 9,
      
      fluidRow(
        
        column(
          width = 7,
          h4("Scatter plot — click an event"),
          plotlyOutput(
            "scatter",
            height = "650px"
          )
        ),
        
        column(
          width = 5,
          h4("Selected event"),
          plotlyOutput(
            "event_trace",
            height = "500px"
          ),
          tableOutput("event_information")
        )
        
      )
      
    )
    
  )
  
)


# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {
  
  selected_event <- reactiveVal(events$event_id[1])
  
  
  # ==========================================================
  # FILTER EVENTS
  # use real data limits with a tiny tolerance
  # ==========================================================
  
  filtered_events <- reactive({
    
    tol <- 1e-9
    
    events |>
      filter(
        event_I_over_I0_pct >= input$ratio_filter[1] - tol,
        event_I_over_I0_pct <= input$ratio_filter[2] + tol,
        dwell_total_ms >= input$dwell_filter[1] - tol,
        dwell_total_ms <= input$dwell_filter[2] + tol
      ) |>
      arrange(event_id)
    
  })
  
  
  # ==========================================================
  # EVENT COUNT
  # ==========================================================
  
  output$event_count <- renderText({
    
    df <- filtered_events()
    
    paste0(
      "Events displayed: ",
      nrow(df),
      " / ",
      nrow(events)
    )
    
  })
  
  
  # ==========================================================
  # KEEP SELECTED EVENT VALID AFTER FILTERING
  # ==========================================================
  
  observe({
    
    df <- filtered_events()
    
    if (nrow(df) == 0) {
      return()
    }
    
    if (!(selected_event() %in% df$event_id)) {
      selected_event(df$event_id[1])
    }
    
  })
  
  
  # ==========================================================
  # CLICK ON SCATTER POINT
  # ==========================================================
  
  observeEvent(
    event_data("plotly_click", source = "event_scatter"),
    {
      
      click <- event_data("plotly_click", source = "event_scatter")
      
      if (!is.null(click$key)) {
        selected_event(as.integer(click$key))
      }
      
    }
  )
  
  
  # ==========================================================
  # GO DIRECTLY TO EVENT
  # ==========================================================
  
  observeEvent(
    input$go_event,
    {
      
      target <- as.integer(input$event_number)
      valid_ids <- filtered_events()$event_id
      
      if (target %in% valid_ids) {
        selected_event(target)
      }
      
    }
  )
  
  
  # ==========================================================
  # PREVIOUS EVENT
  # ==========================================================
  
  observeEvent(
    input$previous_event,
    {
      
      ids <- filtered_events()$event_id
      
      if (length(ids) == 0) {
        return()
      }
      
      current_index <- match(selected_event(), ids)
      
      if (is.na(current_index)) {
        selected_event(ids[1])
      } else if (current_index > 1) {
        selected_event(ids[current_index - 1])
      }
      
    }
  )
  
  
  # ==========================================================
  # NEXT EVENT
  # ==========================================================
  
  observeEvent(
    input$next_event,
    {
      
      ids <- filtered_events()$event_id
      
      if (length(ids) == 0) {
        return()
      }
      
      current_index <- match(selected_event(), ids)
      
      if (is.na(current_index)) {
        selected_event(ids[1])
      } else if (current_index < length(ids)) {
        selected_event(ids[current_index + 1])
      }
      
    }
  )
  
  
  # ==========================================================
  # SCATTER PLOT
  # ==========================================================
  
  output$scatter <- renderPlotly({
    
    df <- filtered_events()
    
    validate(
      need(
        nrow(df) > 0,
        "No events match these filters."
      )
    )
    
    p <- plot_ly(
      data = df,
      x = ~event_I_over_I0_pct,
      y = ~dwell_total_ms,
      key = ~event_id,
      source = "event_scatter",
      type = "scattergl",
      mode = "markers",
      text = ~paste0(
        "<b>Event ", event_id, "</b>",
        "<br>I/I0: ", round(event_I_over_I0_pct, 1), "%",
        "<br>Dwell: ", round(dwell_total_ms, 3), " ms",
        "<br>Baseline: ", round(baseline_pA, 2), " pA"
      ),
      hoverinfo = "text",
      marker = list(
        size = 7,
        opacity = 0.7,
        color = "black"
      )
    )
    
    selected_df <- df |>
      filter(event_id == selected_event())
    
    if (nrow(selected_df) == 1) {
      p <- p |>
        add_trace(
          data = selected_df,
          x = ~event_I_over_I0_pct,
          y = ~dwell_total_ms,
          type = "scatter",
          mode = "markers",
          inherit = FALSE,
          showlegend = FALSE,
          hoverinfo = "skip",
          marker = list(
            size = 15,
            color = "black",
            symbol = "circle-open",
            line = list(
              color = "black",
              width = 3
            )
          )
        )
    }
    
    p |>
      layout(
        showlegend = FALSE,
        xaxis = list(
          title = "I/I0 (%)",
          range = c(
            max(0, RATIO_MIN_UI - 2),
            min(110, RATIO_MAX_UI + 2)
          ),
          zeroline = FALSE
        ),
        yaxis = list(
          title = "Dwell time (ms)",
          type = "log",
          tickmode = "array",
          tickvals = dwell_ticks,
          ticktext = dwell_tick_text,
          zeroline = FALSE
        ),
        hovermode = "closest"
      )
    
  })
  
  
  # ==========================================================
  # SELECTED TRACE DATA
  # ==========================================================
  
  selected_trace_data <- reactive({
    
    traces |>
      filter(event_id == selected_event()) |>
      arrange(time_relative_ms)
    
  })
  
  
  selected_event_data <- reactive({
    
    events |>
      filter(event_id == selected_event()) |>
      slice(1)
    
  })
  
  
  # ==========================================================
  # EVENT TRACE
  # ==========================================================
  
  output$event_trace <- renderPlotly({
    
    tr <- selected_trace_data()
    ev <- selected_event_data()
    
    validate(
      need(
        nrow(tr) > 0,
        "Select an event."
      )
    )
    
    if (input$trace_mode == "ratio") {
      
      y_values <- tr$I_over_I0_trace_pct
      y_title <- "I/I0 (%)"
      baseline_value <- 100
      
    } else {
      
      y_values <- tr$current_pA
      y_title <- "Current (pA)"
      baseline_value <- ev$baseline_pA
      
    }
    
    y_min <- min(y_values, na.rm = TRUE)
    y_max <- max(y_values, na.rm = TRUE)
    
    padding_y <- max(
      1,
      0.05 * (y_max - y_min)
    )
    
    p <- plot_ly(
      x = tr$time_relative_ms,
      y = y_values,
      type = "scatter",
      mode = "lines",
      line = list(
        width = 1.5,
        color = "black"
      ),
      text = paste0(
        "Time: ", round(tr$time_relative_ms, 3), " ms",
        "<br>",
        if (input$trace_mode == "ratio") {
          paste0("I/I0: ", round(y_values, 2), "%")
        } else {
          paste0("Current: ", round(y_values, 2), " pA")
        }
      ),
      hoverinfo = "text"
    )
    
    shapes <- list(
      
      list(
        type = "line",
        x0 = 0,
        x1 = 0,
        y0 = y_min - padding_y,
        y1 = y_max + padding_y,
        line = list(
          color = "grey40",
          dash = "dash",
          width = 1
        )
      ),
      
      list(
        type = "line",
        x0 = ev$dwell_total_ms,
        x1 = ev$dwell_total_ms,
        y0 = y_min - padding_y,
        y1 = y_max + padding_y,
        line = list(
          color = "grey40",
          dash = "dash",
          width = 1
        )
      ),
      
      list(
        type = "line",
        x0 = min(tr$time_relative_ms, na.rm = TRUE),
        x1 = max(tr$time_relative_ms, na.rm = TRUE),
        y0 = baseline_value,
        y1 = baseline_value,
        line = list(
          color = "grey60",
          dash = "dot",
          width = 1
        )
      )
      
    )
    
    p |>
      layout(
        title = list(
          text = paste0(
            "Event ", ev$event_id,
            "<br>",
            round(ev$dwell_total_ms, 3),
            " ms | I/I0 = ",
            round(ev$event_I_over_I0_pct, 1),
            "%"
          )
        ),
        xaxis = list(
          title = "Time relative to event start (ms)",
          zeroline = FALSE
        ),
        yaxis = list(
          title = y_title,
          range = c(
            y_min - padding_y,
            y_max + padding_y
          )
        ),
        shapes = shapes,
        showlegend = FALSE
      )
    
  })
  
  
  # ==========================================================
  # EVENT INFORMATION
  # ==========================================================
  
  output$event_information <- renderTable({
    
    ev <- selected_event_data()
    
    data.frame(
      Parameter = c(
        "Event ID",
        "Dwell time",
        "I/I0",
        "Baseline",
        "Channel",
        "Sweep"
      ),
      Value = c(
        ev$event_id,
        paste0(round(ev$dwell_total_ms, 3), " ms"),
        paste0(round(ev$event_I_over_I0_pct, 1), "%"),
        paste0(round(ev$baseline_pA, 2), " pA"),
        ev$channel,
        ev$sweep
      ),
      check.names = FALSE
    )
    
  })
  
}


# ============================================================
# LAUNCH APP
# ============================================================

shinyApp(
  ui = ui,
  server = server
)