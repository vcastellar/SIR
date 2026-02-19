app_server <- function(input, output, session, models) {

  ## ---------------------------------------------------------
  ## Valores reactivos auxiliares
  ## ---------------------------------------------------------
  rv <- shiny::reactiveValues(
    derived_selected = NULL
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
      choices  = model()$states,
      selected = model()$states,
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
      model()$states,
      function(s) {
        id <- paste0("init_", s)
        shiny::req(!is.null(input[[id]]))
        input[[id]]
      },
      numeric(1)
    )

    shiny::req(!any(is.na(vals)))
    setNames(as.numeric(vals), model()$states)
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
  ## Gestión de selección de variables derivadas (sticky)
  ## ---------------------------------------------------------
  observeEvent(sim(), {
    derived <- sim()$derived
    if (is.null(derived)) derived <- sim()$flows
    if (is.null(derived) || ncol(derived) <= 1) return()

    derived_names <- setdiff(names(derived), "time")

    if (is.null(rv$derived_selected)) {
      # primera vez
      rv$derived_selected <- derived_names
    } else {
      # mantener solo los que siguen existiendo
      rv$derived_selected <- intersect(rv$derived_selected, derived_names)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$derived_to_plot, {
    rv$derived_selected <- input$derived_to_plot
  }, ignoreInit = TRUE)

  ## ---------------------------------------------------------
  ## Selector de variables derivadas
  ## ---------------------------------------------------------
  output$derived_select_ui <- shiny::renderUI({
    shiny::req(sim())

    derived <- sim()$derived
    if (is.null(derived)) derived <- sim()$flows
    if (is.null(derived) || ncol(derived) <= 1) {
      return(shiny::tags$em("No derived variables defined for this model"))
    }

    derived_names <- setdiff(names(derived), "time")

    shiny::checkboxGroupInput(
      "derived_to_plot",
      NULL,
      choices  = derived_names,
      selected = rv$derived_selected,
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
    sim2$model$states <- states

    plot(sim2, what = "states")
  })

  ## ---------------------------------------------------------
  ## Plot de variables derivadas
  ## ---------------------------------------------------------
  output$plot_derived <- shiny::renderPlot({
    shiny::req(sim(), input$derived_to_plot)

    derived <- input$derived_to_plot
    if (length(derived) == 0) {
      plot.new()
      text(0.5, 0.5, "No derived variables selected")
      return()
    }

    sim2 <- sim()
    keep <- c("time", derived)

    derived_data <- sim2$derived
    if (is.null(derived_data)) derived_data <- sim2$flows

    sim2$states <- derived_data[, keep, drop = FALSE]
    sim2$model$states <- derived

    plot(sim2, what = "states")
  })
}

