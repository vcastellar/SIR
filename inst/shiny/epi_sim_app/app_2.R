## ============================================================================
## Epidemic model simulator (ODE-based, Shiny app)
## ============================================================================

library(shiny)
library(SIR)

## ---------------------------------------------------------------------------
## Utilidad: operador NULL-coalescing
## ---------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

## ---------------------------------------------------------------------------
## Helper: sliders de parámetros
## ---------------------------------------------------------------------------
param_sliders_ui <- function(model) {

  lapply(model$par_names, function(p) {

    val   <- model$defaults[[p]] %||% 1
    lower <- model$lower[[p]]    %||% (val / 10)
    upper <- model$upper[[p]]    %||% (val * 10)

    sliderInput(
      inputId = paste0("par_", p),
      label   = p,
      min     = lower,
      max     = upper,
      value   = val
    )
  })
}

## ---------------------------------------------------------------------------
## Helper: sliders de estados iniciales
## ---------------------------------------------------------------------------
init_sliders_ui <- function(model) {

  lapply(model$state_names, function(s) {

    val <- model$init[[s]] %||% 0

    sliderInput(
      inputId = paste0("init_", s),
      label   = s,
      min     = 0,
      max     = max(1, val * 5),
      value   = val
    )
  })
}

## ---------------------------------------------------------------------------
## Helper: mostrar RHS
## ---------------------------------------------------------------------------
rhs_text <- function(model) {
  paste(deparse(body(model$rhs)), collapse = "\n")
}

## ---------------------------------------------------------------------------
## UI factory
## ---------------------------------------------------------------------------
app_ui <- function(models) {

  fluidPage(

    titlePanel("Epidemic model simulator (ODE-based)"),

    sidebarLayout(

      sidebarPanel(

        selectInput(
          "model_name",
          "Epidemic model",
          choices  = names(models),
          selected = names(models)[1]
        ),

        hr(),
        h4("Parameters"),
        uiOutput("params_ui"),

        hr(),
        h4("Initial conditions"),
        uiOutput("init_ui"),

        hr(),
        h4("States to plot"),
        uiOutput("states_select_ui"),

        hr(),
        h4("Flows to plot"),
        uiOutput("flows_select_ui"),

        hr(),
        sliderInput(
          "t_max",
          "Time horizon",
          min   = 10,
          max   = 500,
          value = 150
        )
      ),

      mainPanel(
        tabsetPanel(
          tabPanel(
            "States",
            plotOutput("plot_states", height = 400)
          ),
          tabPanel(
            "Flows",
            plotOutput("plot_flows", height = 400)
          ),
          tabPanel(
            "Equations",
            tags$pre(
              style = "font-size: 13px;",
              verbatimTextOutput("rhs_equations")
            )
          )
        )
      )
    )
  )
}

## ---------------------------------------------------------------------------
## Server factory
## ---------------------------------------------------------------------------
app_server <- function(input, output, session, models) {

  ## ---------------------------------------------------------
  ## Modelo seleccionado
  ## ---------------------------------------------------------
  model <- reactive({
    req(input$model_name)
    models[[input$model_name]]
  })

  ## ---------------------------------------------------------
  ## UI dinámica
  ## ---------------------------------------------------------
  output$params_ui <- renderUI({
    tagList(param_sliders_ui(model()))
  })

  output$init_ui <- renderUI({
    tagList(init_sliders_ui(model()))
  })

  output$states_select_ui <- renderUI({
    checkboxGroupInput(
      "states_to_plot",
      NULL,
      choices  = model()$state_names,
      selected = model()$state_names,
      inline   = TRUE
    )
  })

  output$rhs_equations <- renderText({
    rhs_text(model())
  })

  ## ---------------------------------------------------------
  ## Parámetros
  ## ---------------------------------------------------------
  parms <- reactive({
    req(model())

    vals <- sapply(
      model()$par_names,
      function(p) input[[paste0("par_", p)]]
    )

    req(!any(is.null(vals)))
    vals
  })

  ## ---------------------------------------------------------
  ## Estados iniciales
  ## ---------------------------------------------------------
  init <- reactive({
    req(model())

    vals <- sapply(
      model()$state_names,
      function(s) input[[paste0("init_", s)]]
    )

    req(!any(is.null(vals)))
    vals
  })

  ## ---------------------------------------------------------
  ## Simulación
  ## ---------------------------------------------------------
  sim <- reactive({
    req(parms(), init(), input$t_max)

    simulate_epi(
      model = model(),
      times = 0:input$t_max,
      parms = parms(),
      init  = init()
    )
  })

  ## ---------------------------------------------------------
  ## Selector de flujos
  ## ---------------------------------------------------------
  output$flows_select_ui <- renderUI({
    req(sim())

    flows <- sim()$flows

    if (is.null(flows) || ncol(flows) <= 1) {
      return(tags$em("No flows defined for this model"))
    }

    flow_names <- setdiff(names(flows), "time")

    checkboxGroupInput(
      "flows_to_plot",
      NULL,
      choices  = flow_names,
      selected = flow_names,
      inline   = TRUE
    )
  })

  ## ---------------------------------------------------------
  ## Plot de estados
  ## ---------------------------------------------------------
  output$plot_states <- renderPlot({
    req(sim(), input$states_to_plot)

    states <- input$states_to_plot

    if (length(states) == 0) {
      plot.new()
      text(0.5, 0.5, "No states selected")
      return()
    }

    sim2 <- sim()
    keep <- c("time", states)

    sim2$states <- sim2$states[, keep, drop = FALSE]
    sim2$model$state_names <- states

    plot(sim2, what = "states")
  })

  ## ---------------------------------------------------------
  ## Plot de flujos
  ## ---------------------------------------------------------
  output$plot_flows <- renderPlot({
    req(sim(), input$flows_to_plot)

    flows <- input$flows_to_plot

    if (length(flows) == 0) {
      plot.new()
      text(0.5, 0.5, "No flows selected")
      return()
    }

    sim2 <- sim()
    keep <- c("time", flows)

    sim2$states <- sim2$flows[, keep, drop = FALSE]
    sim2$model$state_names <- flows

    plot(sim2, what = "states")
  })
}
