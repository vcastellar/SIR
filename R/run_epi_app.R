#' Run the epidemic model simulator
#'
#' Launches the interactive Shiny application for exploring and simulating
#' registered epidemiological models.
#'
#' The app automatically detects all models available through
#' \code{list_epi_models()}, including user-defined models registered with
#' \code{register_epi_model()}.
#'
#' @export
run_epi_app <- function() {


  app_dir <- system.file(
    "shiny", "epi_sim_app",
    package = "SIR"
  )

  if (app_dir == "") {
    stop("Could not find Shiny app directory.", call. = FALSE)
  }

  shiny::runApp(app_dir)
}


