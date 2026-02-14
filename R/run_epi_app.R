#' Launch the epidemic model Shiny simulator
#'
#' Launches the interactive Shiny application for simulating deterministic
#' epidemiological models defined as systems of ordinary differential equations
#' (ODEs).
#'
#' By default, the app includes a set of built-in models (e.g. SI, SIR, SEIR).
#' User-defined models can be explored by passing them explicitly via the
#' \code{models} argument.
#'
#' @param models Optional named list of objects of class \code{epi_model}.
#'   These models are added to the built-in ones and become available in the
#'   Shiny interface. If \code{NULL} (default), only built-in models are used.
#'
#' @details
#' The Shiny application does not rely on global state or registries.
#' All models must be provided explicitly at launch time so that they are
#' available in the Shiny session.
#'
#' User-defined models must be fully specified \code{epi_model} objects,
#' including state variables and parameter definitions.
#'
#' @return
#' A \code{shiny.appobj}. This function is called for its side effects
#' (launching the Shiny application).
#'
#' @examples
#' ## ---------------------------------------------------------
#' ## Define a custom SEIRD model
#' ## ---------------------------------------------------------
#'
#' seird_rhs <- function(time, state, parms) {
#'   with(as.list(c(state, parms)), {
#'
#'     N <- S + E + I + R
#'
#'     lambda <- beta * S * I / N
#'
#'     dS <- -lambda
#'     dE <-  lambda - sigma * E
#'     dI <-  sigma * E - gamma * I - mu * I
#'     dR <-  gamma * I
#'     dD <-  mu * I
#'
#'     list(
#'       c(dS, dE, dI, dR, dD),
#'       incidence = lambda,
#'       deaths    = mu * I
#'     )
#'   })
#' }
#'
#' seird_model <- epi_model(
#'   name        = "SEIRD",
#'   rhs         = seird_rhs,
#'   states = c("S", "E", "I", "R", "D"),
#'   flows       = c("incidence", "deaths"),
#'   par_names   = c("beta", "sigma", "gamma", "mu"),
#'   defaults    = c(beta = 0.4, sigma = 0.2, gamma = 0.1, mu = 0.02)
#' )
#'
#' ## ---------------------------------------------------------
#' ## Launch the Shiny app with the custom model
#' ## ---------------------------------------------------------
#'
#' run_epi_app(
#'   models = list(SEIRD = seird_model)
#' )
#'
#' @seealso
#' \code{\link{epi_model}}, \code{\link{simulate_epi}}
#'
#' @export
run_epi_app <- function(models = NULL) {

  # builtin <- .get_builtin_models()
  builtin <- setNames(
    lapply(list_models(), get_model),
    list_models()
  )

  if (is.null(models)) {
    all_models <- builtin
  } else {
    stopifnot(is.list(models))
    stopifnot(length(models) > 0)

    ok <- vapply(models, inherits, logical(1), "epi_model")
    if (!all(ok)) {
      stop("All elements of 'models' must be objects of class 'epi_model'.")
    }

    if (is.null(names(models)) || any(names(models) == "")) {
      stop("'models' must be a named list.")
    }

    dup <- intersect(names(models), names(builtin))
    if (length(dup) > 0) {
      stop(
        "Model name(s) already exist: ",
        paste(dup, collapse = ", ")
      )
    }

    all_models <- c(builtin, models)
  }

  shiny::shinyApp(
    ui = app_ui(all_models),
    server = function(input, output, session) {
      app_server(input, output, session, all_models)
    }
  )
}

