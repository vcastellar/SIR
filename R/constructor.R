#-------------------------------------------------------------------------------
# Epidemiological model constructor
#-------------------------------------------------------------------------------
#' Construct an epidemiological model
#'
#' Creates an \code{epi_model} object describing the structure of a deterministic
#' compartmental epidemiological model. The model definition includes the
#' system of equations, state variables, optional derived variables, parameter
#' names.
#'
#' The resulting object is a *model definition* and does not perform any
#' simulation by itself. It is intended to be used by downstream functions such
#' as \code{simulate_epi()}.
#'
#' @param name Character scalar giving the name of the model.
#'
#' @param rhs A function defining the right-hand side of the system of
#'   differential (or difference) equations. The function must be compatible
#'   with the simulation backend used by \code{simulate_epi()}.
#'
#' @param par_names Character vector with the names of the model parameters.
#'
#' @param states Character vector giving the names of the compartmental state
#'   variables of the model. Each state name must be unique.
#'
#' @param flows Optional character vector giving the names of flow variables
#'   used internally by the model or returned by the RHS. Deprecated in favor
#'   of \code{derived}; kept for backward compatibility.
#'
#' @param derived Optional character vector giving the names of derived
#'   variables returned by the RHS and extracted by \code{simulate_epi()}.
#'
#' @param lower Optional named numeric vector giving lower bounds for model
#'   parameters. Names must match \code{par_names}.
#'
#' @param upper Optional named numeric vector giving upper bounds for model
#'   parameters. Names must match \code{par_names}.
#'
#' @param defaults Optional named numeric vector of default parameter values.
#'   Names must match \code{par_names}.
#'
#' @param init Optional named numeric vector of initial values for the state
#'   variables. Names must match \code{states}.
#'
#' @return
#' An object of class \code{"epi_model"} containing the model definition.
#'
#' @seealso
#' \code{\link{simulate_epi}}
#'
#' @export

epi_model <- function(name,
                      rhs,
                      par_names,
                      states,
                      flows = character(0),
                      derived = NULL,
                      lower = NULL,
                      upper = NULL,
                      defaults = NULL,
                      init = NULL) {

  ## --- basic checks ----------------------------------------------------------
  stopifnot(is.character(name), length(name) == 1)
  stopifnot(is.function(rhs))
  stopifnot(is.character(par_names), length(par_names) >= 1)
  stopifnot(is.character(states), length(states) >= 1)
  stopifnot(length(unique(states)) == length(states))

  using_derived <- !is.null(derived)
  if (using_derived) {
    if (!missing(flows) && length(flows) > 0) {
      stop("Use either `flows` (deprecated) or `derived`, not both.")
    }
    flows <- derived
  }

  if (missing(flows) || length(flows) == 0) {
    flows <- character(0)
  } else if (!using_derived) {
    .Deprecated(msg = paste0(
      "`flows` in `epi_model()` is deprecated and will be removed in a future release. ",
      "Use `derived` instead."
    ))
  }

  ## a variable cannot be both state and flow
  stopifnot(
    is.character(flows),
    length(unique(flows)) == length(flows),
    !any(flows %in% states)
  )

  #-----------------------------------------------------------------------------
  ## Parameter bounds
  #-----------------------------------------------------------------------------
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

  ## --- defaults --------------------------------------------------------------
  if (!is.null(defaults)) {
    stopifnot(is.numeric(defaults), all(par_names %in% names(defaults)))
    defaults <- defaults[par_names]
  }

  ## --- init ------------------------------------------------------------------
  if (!is.null(init)) {
    stopifnot(is.numeric(init), all(states %in% names(init)))
    init <- init[states]
  }

  structure(
    list(
      name = name,
      rhs = rhs,
      par_names = par_names,
      states = states,
      derived = flows,
      flows = flows,
      lower = lower,
      upper = upper,
      defaults = defaults,
      init = init
    ),
    class = "epi_model"
  )
}


#' Print an epidemic model object
#'
#' @name print.epi_model
#' @description
#' Provides a concise, human-readable summary of an \code{epi_model} object.
#' The printed output includes the model name, state variables, parameters,
#' the declared model derived variables, and the underlying system of differential
#' equations.
#'
#' This method is automatically called when an object of class
#' \code{"epi_model"} is printed at the console.
#'
#' @details
#' Model derived variables are shown when available.
#' These may include incidence, recoveries, or any other observable
#' returned by the model's right-hand side (\code{rhs}) function.
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

  cat("<epi_model> ", x$name, "\n", sep = "")

  ## --- core structure --------------------------------------------------------
  cat("  States:   ", paste(x$states, collapse = ", "), "\n", sep = "")
  cat("  Params:   ", paste(x$par_names, collapse = ", "), "\n", sep = "")

  ## --- derived variables -----------------------------------------------------
  derived_names <- x$derived
  if (is.null(derived_names)) derived_names <- x$flows

  if (!is.null(derived_names) && length(derived_names) > 0) {
    cat("  Derived:  ", paste(derived_names, collapse = ", "), "\n", sep = "")
  }

  ## --- parameter bounds ------------------------------------------------------
  if (!is.null(x$lower) && !is.null(x$upper)) {
    cat("  Bounds:\n")
    b <- cbind(lower = x$lower, upper = x$upper)
    print(b)
  }

  ## --- equations -------------------------------------------------------------
  if (!is.null(x$rhs)) {
    cat("  Equations (rhs):\n")
    rhs_txt <- deparse(x$rhs)
    rhs_txt <- rhs_txt[nzchar(trimws(rhs_txt))]
    for (ln in rhs_txt) {
      cat("    ", ln, "\n", sep = "")
    }
  }

  invisible(x)
}
