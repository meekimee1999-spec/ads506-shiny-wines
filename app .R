library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(tsibble)
library(fable)
library(fabletools)
library(feasts)
library(ggplot2)
library(urca)
library(stringr)

# Load and wrangle AustralianWines.csv into a monthly tsibble
raw_wines <- read_csv("AustralianWines.csv", show_col_types = FALSE)

wine_ts <- raw_wines |>
  mutate(Month = yearmonth(parse_date_time(Month, orders = c("my", "ym", "mdy", "ymd")))) |>
  pivot_longer(-Month, names_to = "Varietal", values_to = "Sales") |>
  mutate(Sales = as.numeric(Sales)) |>
  filter(!is.na(Month), !is.na(Sales)) |>
  as_tsibble(index = Month, key = Varietal)

varietals_all <- sort(unique(wine_ts$Varietal))
min_date <- as.Date(min(wine_ts$Month))
max_date <- as.Date(max(wine_ts$Month))

ui <- page_sidebar(
  title = "Australian Wine Sales — Forecasting Explorer",
  theme = bs_theme(bootswatch = "flatly"),
  sidebar = sidebar(
    width = 300,
    selectizeInput(
      "varietals",
      "Varietals (select 1 or more)",
      choices = varietals_all,
      selected = head(varietals_all, 2),
      multiple = TRUE
    ),
    dateRangeInput(
      "daterange",
      "Date range",
      start = min_date,
      end = max_date,
      min = min_date,
      max = max_date
    ),
    sliderInput("val_months", "Validation window (months)", min = 6, max = 36, value = 12, step = 6),
    sliderInput("horizon", "Forecast horizon (months)", min = 6, max = 36, value = 12, step = 6),
    sliderInput("hist_years", "History shown on forecast plot (years)", min = 2, max = 15, value = 5, step = 1),
    checkboxInput("log_transform", "Log-transform Sales (stabilize seasonal variance)", value = FALSE),
    checkboxInput("use_ensemble", "Include Ensemble (TSLM + ETS + ARIMA)", value = TRUE),
    actionButton("fit_btn", "Fit models", class = "btn-primary w-100")
  ),
  navset_card_tab(
    nav_panel(
      "Overview",
      plotOutput("overview_plot", height = "340px"),
      plotOutput("stl_plot", height = "380px")
    ),
    nav_panel(
      "Forecasts",
      plotOutput("fc_plot", height = "520px")
    ),
    nav_panel(
      "Model specs",
      tableOutput("specs_table")
    ),
    nav_panel(
      "Accuracy",
      h5("Validation Accuracy (Holdout Window)"),
      tableOutput("val_acc_table"),
      hr(),
      h5("Training Accuracy"),
      tableOutput("train_acc_table")
    ),
    nav_panel(
      "About",
      h4("Data Source & Modeling Choices"),
      p("Dataset: AustralianWines.csv (monthly Australian wine sales by varietal)."),
      p("Models: TSLM (trend + season), automatic ETS, and automatic ARIMA fit per varietal on the training split (all filtered data prior to the validation holdout window)."),
      p("Reproducibility: Select varietals, date range, validation window, and horizon in the sidebar, click 'Fit models', and inspect the Forecasts, Model specs, and Accuracy tabs."),
      p(em("AI disclosure: App structure drafted with AI assistance and verified/deployed in RStudio."))
    )
  )
)

server <- function(input, output, session) {
  filtered_data <- reactive({
    req(input$varietals, input$daterange)
    wine_ts |>
      filter(
        Varietal %in% input$varietals,
        Month >= yearmonth(input$daterange[1]),
        Month <= yearmonth(input$daterange[2])
      )
  })
  
  output$overview_plot <- renderPlot({
    df <- filtered_data()
    req(nrow(df) > 0)
    ggplot(df, aes(x = Month, y = Sales, color = Varietal)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~Varietal, scales = "free_y", ncol = 1) +
      labs(title = "Monthly Sales by Varietal", x = NULL, y = "Sales") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")
  })
  
  output$stl_plot <- renderPlot({
    df <- filtered_data()
    req(nrow(df) > 0)
    first_var <- input$varietals[1]
    df |>
      filter(Varietal == first_var) |>
      model(STL(Sales ~ season(window = "periodic"))) |>
      components() |>
      autoplot() +
      labs(title = paste("STL Decomposition —", first_var)) +
      theme_minimal(base_size = 11)
  })
  
  fitted_results <- eventReactive(input$fit_btn, {
    df <- filtered_data()
    req(nrow(df) > 0)
    
    max_m <- max(df$Month)
    train_end <- max_m - input$val_months
    train_df <- df |> filter(Month <= train_end)
    val_df <- df |> filter(Month > train_end)
    
    if (isTRUE(input$log_transform)) {
      fits <- train_df |>
        model(
          TSLM = TSLM(log(Sales) ~ trend() + season()),
          ETS = ETS(log(Sales)),
          ARIMA = ARIMA(log(Sales))
        )
    } else {
      fits <- train_df |>
        model(
          TSLM = TSLM(Sales ~ trend() + season()),
          ETS = ETS(Sales),
          ARIMA = ARIMA(Sales)
        )
    }
    
    if (isTRUE(input$use_ensemble)) {
      fits <- fits |> mutate(Ensemble = (TSLM + ETS + ARIMA) / 3)
    }
    
    val_fc <- fits |> forecast(new_data = val_df)
    fut_fc <- fits |> forecast(h = input$horizon)
    
    list(
      df = df,
      train_df = train_df,
      val_df = val_df,
      fits = fits,
      val_fc = val_fc,
      fut_fc = fut_fc
    )
  }, ignoreNULL = FALSE)
  
  output$fc_plot <- renderPlot({
    res <- fitted_results()
    hist_start <- max(res$df$Month) - (input$hist_years * 12)
    hist_df <- res$df |> filter(Month >= hist_start)
    
    res$fut_fc |>
      autoplot(hist_df, level = 95, alpha = 0.25) +
      facet_wrap(~Varietal, scales = "free_y", ncol = 1) +
      labs(
        title = "Forecast Comparison with 95% Prediction Intervals",
        x = NULL,
        y = "Sales"
      ) +
      theme_minimal(base_size = 12)
  })
  
  output$specs_table <- renderTable({
    res <- fitted_results()
    res$fits |>
      select(Varietal, TSLM, ETS, ARIMA) |>
      mutate(
        TSLM = format(TSLM),
        ETS = format(ETS),
        ARIMA = format(ARIMA)
      )
  })
  
  output$val_acc_table <- renderTable({
    res <- fitted_results()
    res$val_fc |>
      accuracy(res$val_df) |>
      select(Varietal, .model, RMSE, MAE, MAPE) |>
      arrange(Varietal, RMSE)
  })
  
  output$train_acc_table <- renderTable({
    res <- fitted_results()
    res$fits |>
      accuracy() |>
      select(Varietal, .model, RMSE, MAE, MAPE) |>
      arrange(Varietal, RMSE)
  })
}

shinyApp(ui, server)