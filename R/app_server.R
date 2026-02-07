app_server <- function(input, output, session, models) {

  ## ---------------------------------------------------------
  ## Valores reactivos auxiliares
  ## ---------------------------------------------------------
  rv <- shiny::reactiveValues(
    flows_selected = NULL
  )

  ## ---------------------------------------------------------
  ## Modelo seleccionado
  ## ---------------------------------------------------------
  model <- shiny::reactive({
    shiny::req(input$model_name)
    models[[input$model_name]]
  })

  ## ---------------------------------------------------------
  ## UI dinámica
  ## ---------------------------------------------------------
  output$params_ui <- shiny::renderUI({
    shiny::tagList(param_sliders_ui(model()))
  })

  output$init_ui <- shiny::renderUI({
    shiny::tagList(init_sliders_ui(model()))
  })

  output$states_select_ui <- shiny::renderUI({
    shiny::checkboxGroupInput(
      "states_to_plot",
      NULL,
      choices  = model()$state_names,
      selected = model()$state_names,
      inline   = TRUE
    )
  })

  output$rhs_equations <- shiny::renderText({
    rhs_text(model())
  })

  ## ---------------------------------------------------------
  ## Parámetros
  ## ---------------------------------------------------------
  parms <- shiny::reactive({
    shiny::req(model())

    vals <- vapply(
      model()$par_names,
      function(p) {
        id <- paste0("par_", p)
        shiny::req(!is.null(input[[id]]))
        input[[id]]
      },
      numeric(1)
    )

    shiny::req(!any(is.na(vals)))
    setNames(as.numeric(vals), model()$par_names)
  })



  ## ---------------------------------------------------------
  ## Estados iniciales
  ## ---------------------------------------------------------
  init <- shiny::reactive({
    shiny::req(model())

    vals <- vapply(
      model()$state_names,
      function(s) {
        id <- paste0("init_", s)
        shiny::req(!is.null(input[[id]]))
        input[[id]]
      },
      numeric(1)
    )

    shiny::req(!any(is.na(vals)))
    setNames(as.numeric(vals), model()$state_names)
  })

  ## ---------------------------------------------------------
  ## Simulación
  ## ---------------------------------------------------------
  sim <- shiny::reactive({
    shiny::req(parms(), init(), input$t_max)

    simulate_epi(
      model = model(),
      times = 0:input$t_max,
      parms = parms(),
      init  = init()
    )
  })

  ## ---------------------------------------------------------
  ## Gestión de selección de flujos (sticky)
  ## ---------------------------------------------------------
  observeEvent(sim(), {
    flows <- sim()$flows
    if (is.null(flows) || ncol(flows) <= 1) return()

    flow_names <- setdiff(names(flows), "time")

    if (is.null(rv$flows_selected)) {
      # primera vez
      rv$flows_selected <- flow_names
    } else {
      # mantener solo los que siguen existiendo
      rv$flows_selected <- intersect(rv$flows_selected, flow_names)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$flows_to_plot, {
    rv$flows_selected <- input$flows_to_plot
  }, ignoreInit = TRUE)

  ## ---------------------------------------------------------
  ## Selector de flujos
  ## ---------------------------------------------------------
  output$flows_select_ui <- shiny::renderUI({
    shiny::req(sim())

    flows <- sim()$flows
    if (is.null(flows) || ncol(flows) <= 1) {
      return(shiny::tags$em("No flows defined for this model"))
    }

    flow_names <- setdiff(names(flows), "time")

    shiny::checkboxGroupInput(
      "flows_to_plot",
      NULL,
      choices  = flow_names,
      selected = rv$flows_selected,
      inline   = TRUE
    )
  })

  ## ---------------------------------------------------------
  ## Plot de estados
  ## ---------------------------------------------------------
  output$plot_states <- shiny::renderPlot({
    shiny::req(sim(), input$states_to_plot)

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
  output$plot_flows <- shiny::renderPlot({
    shiny::req(sim(), input$flows_to_plot)

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

