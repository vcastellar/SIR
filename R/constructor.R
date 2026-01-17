#' Create a new epidemic model object
#'
#' @name new_epi_model
#' @description
#' Constructs an object of class \code{"epi_model"} representing a deterministic
#' compartmental epidemic model. In addition to the ODE system, the model must
#' explicitly declare which output of the right-hand side corresponds to the
#' latent incidence rate.
#'
#' @details
#' ## Incidence
#' The \code{rhs} function may return additional named outputs besides the state
#' derivatives. One of these outputs must represent the **latent incidence rate**
#' (e.g. force of infection or case-generation rate).
#'
#' The \code{incidence} argument specifies the name of this output. This makes the
#' definition of epidemic cases explicit and allows generic functions such as
#' \code{simulate_epi()}, \code{fit_epi_model()}, and \code{summary()} to operate
#' consistently across models.
#'
#' @param name Character scalar. Human-readable name of the model (e.g. \code{"SIR"}).
#' @param rhs Function defining the ODE system. Must have signature
#'   \code{function(time, state, parms)} and return a list whose first element is
#'   the vector of state derivatives.
#' @param state_names Character vector of state variable names.
#' @param par_names Character vector of parameter names.
#' @param incidence Character scalar. Name of the \code{rhs} output corresponding
#'   to the latent incidence rate.
#' @param lower Optional named numeric vector of lower parameter bounds.
#' @param upper Optional named numeric vector of upper parameter bounds.
#' @param defaults Optional named numeric vector of default parameter values.
#' @param init Optional named numeric vector of default initial conditions.
#'
#' @return
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' sir_model <- new_epi_model(
#'   name = "SIR",
#'   rhs = sir_rhs,
#'   state_names = c("S", "I", "R"),
#'   par_names = c("beta", "gamma"),
#'   incidence = "incidence"
#' )
#'
#' @export
new_epi_model <- function(name,
                          rhs,
                          state_names,
                          par_names,
                          incidence,
                          lower = NULL,
                          upper = NULL,
                          defaults = NULL,
                          init = NULL) {

  stopifnot(is.character(name), length(name) == 1)
  stopifnot(is.function(rhs))
  stopifnot(is.character(state_names), length(state_names) >= 1)
  stopifnot(is.character(par_names), length(par_names) >= 1)

  ## --- incidence is mandatory -------------------------------------------------
  if (missing(incidence)) {
    stop("`incidence` must be provided and must name an output returned by `rhs`.")
  }
  stopifnot(is.character(incidence), length(incidence) == 1)

  ## --- parameter bounds -------------------------------------------------------
  if (!is.null(lower)) {
    stopifnot(is.numeric(lower), all(par_names %in% names(lower)))
    lower <- lower[par_names]
  }
  if (!is.null(upper)) {
    stopifnot(is.numeric(upper), all(par_names %in% names(upper)))
    upper <- upper[par_names]
  }
  if (!is.null(lower) && !is.null(upper)) {
    if (any(lower >= upper)) stop("Invalid bounds: lower >= upper.")
  }

  ## --- defaults ---------------------------------------------------------------
  if (!is.null(defaults)) {
    stopifnot(is.numeric(defaults), all(par_names %in% names(defaults)))
    defaults <- defaults[par_names]
  }

  ## --- init -------------------------------------------------------------------
  if (!is.null(init)) {
    stopifnot(is.numeric(init), all(state_names %in% names(init)))
    init <- init[state_names]
  }

  structure(
    list(
      name = name,
      rhs = rhs,
      state_names = state_names,
      par_names = par_names,
      incidence = incidence,
      lower = lower,
      upper = upper,
      defaults = defaults,
      init = init
    ),
    class = "epi_model"
  )
}


#' Print method for epidemic model objects
#' Print an epidemic model object
#'
#' @name print.epi_model
#' @description
#' Provides a concise, human-readable summary of an \code{epi_model} object.
#' The printed output includes the model name, state variables, parameters,
#' the declared definition of epidemic incidence, and the underlying system
#' of differential equations.
#'
#' This method is automatically called when an object of class
#' \code{"epi_model"} is printed at the console.
#'
#' @details
#' The \code{epi_model} class requires an explicit declaration of epidemic
#' incidence via the \code{incidence} field, which specifies the name of the
#' output returned by the model's right-hand side (\code{rhs}) corresponding
#' to the latent incidence rate. This information is displayed by the
#' \code{print()} method to make the epidemiological interpretation of the
#' model explicit.
#'
#' Parameter bounds are shown when available. The model equations are printed
#' by deparsing the \code{rhs} function for inspection.
#'
#' @param x An object of class \code{"epi_model"}.
#' @param ... Further arguments (ignored).
#'
#' @return
#' Invisibly returns the input object \code{x}.
#'
#' @examples
#' SIR_MODEL
#'
#' @export
print.epi_model <- function(x, ...) {

  stopifnot(inherits(x, "epi_model"))

  cat("<epi_model>", x$name, "\n")

  cat("  States:    ", paste(x$state_names, collapse = ", "), "\n", sep = "")
  cat("  Params:    ", paste(x$par_names, collapse = ", "), "\n", sep = "")

  ## --- incidence --------------------------------------------------------------
  if (!is.null(x$incidence)) {
    cat("  Incidence: ", x$incidence, " (rhs output)\n", sep = "")
  }

  ## --- parameter bounds -------------------------------------------------------
  if (!is.null(x$lower) && !is.null(x$upper)) {
    cat("  Bounds:\n")
    b <- cbind(lower = x$lower, upper = x$upper)
    print(b)
  }

  ## --- equations --------------------------------------------------------------
  if (!is.null(x$rhs)) {
    cat("  Equations (rhs):\n")
    rhs_txt <- deparse(x$rhs)
    rhs_txt <- rhs_txt[nzchar(trimws(rhs_txt))]
    for (ln in rhs_txt) {
      cat("   ", ln, "\n", sep = "")
    }
  }

  invisible(x)
}

