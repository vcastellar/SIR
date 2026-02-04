app_ui <- function(models) {

  shiny::fluidPage(

    shiny::titlePanel("Epidemic model simulator (ODE-based)"),

    shiny::sidebarLayout(

      shiny::sidebarPanel(

        shiny::selectInput(
          "model_name",
          "Epidemic model",
          choices  = names(models),
          selected = names(models)[1]
        ),

        shiny::hr(),
        shiny::h4("Parameters"),
        shiny::uiOutput("params_ui"),

        shiny::hr(),
        shiny::h4("Initial conditions"),
        shiny::uiOutput("init_ui"),

        shiny::hr(),
        shiny::h4("States to plot"),
        shiny::uiOutput("states_select_ui"),

        shiny::hr(),
        shiny::h4("Flows to plot"),
        shiny::uiOutput("flows_select_ui"),

        shiny::hr(),
        shiny::sliderInput(
          "t_max",
          "Time horizon",
          min   = 10,
          max   = 500,
          value = 150
        )
      ),

      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel(
            "States",
            shiny::plotOutput("plot_states", height = 400)
          ),
          shiny::tabPanel(
            "Flows",
            shiny::plotOutput("plot_flows", height = 400)
          ),
          shiny::tabPanel(
            "Equations",
            shiny::tags$pre(
              style = "font-size: 13px;",
              shiny::verbatimTextOutput("rhs_equations")
            )
          )
        )
      )
    )
  )
}
