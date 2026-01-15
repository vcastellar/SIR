#' Create a new epidemic model object
#'
#' @name new_epi_model
#'
#' @description
#' Constructs an object of class \code{"epi_model"} that encapsulates all the
#' information required to define and simulate a deterministic compartmental
#' epidemic model (such as SIR, SIRS, or user-defined extensions).
#'
#' An \code{epi_model} object stores the model equations (ODE right-hand side),
#' state variables, parameters, parameter bounds, default values, and optional
#' helpers for building initial conditions and interpreting model output.
#'
#' @details
#' The \code{epi_model} class provides a lightweight abstraction layer separating
#' model definition from simulation, inference, and visualization. Generic tools
#' such as \code{\link{simulate_epi}} operate on this object without requiring
#' model-specific logic.
#'
#' ## Required components
#' At minimum, a model must define:
#' \itemize{
#'   \item a name (\code{name});
#'   \item a right-hand side ODE function compatible with \code{deSolve::ode()}
#'     (\code{rhs});
#'   \item the names of the state variables (\code{state_names});
#'   \item the names of the model parameters (\code{par_names}).
#' }
#'
#' ## Optional components
#' Additional components may be provided to facilitate simulation and fitting:
#' \itemize{
#'   \item lower and upper bounds for parameters (\code{lower}, \code{upper});
#'   \item default parameter values (\code{defaults});
#'   \item a function to construct initial conditions (\code{make_init});
#'   \item a character vector of model equations for display (\code{equations});
#'   \item conventions for output columns such as incidence and cumulative counts
#'     (\code{output}).
#' }
#'
#' Parameter bounds and defaults, when supplied, must be named numeric vectors
#' whose names match \code{par_names}.
#'
#' @param name Character scalar. Human-readable name of the model (e.g. \code{"SIR"},
#'   \code{"SIRS"}).
#' @param rhs Function defining the right-hand side of the ODE system. Must have
#'   signature \code{function(time, state, parms)} and be compatible with
#'   \code{deSolve::ode()}.
#' @param state_names Character vector giving the names of the state variables
#'   required by the model.
#' @param par_names Character vector giving the names of the model parameters.
#' @param lower Optional named numeric vector of lower bounds for parameters.
#' @param upper Optional named numeric vector of upper bounds for parameters.
#' @param defaults Optional named numeric vector of default parameter values.
#' @param make_init Optional function to construct the initial state vector.
#'   Typically takes arguments such as population size and initial conditions
#'   (e.g. \code{N}, \code{I0}, \code{R0}) and returns a named numeric vector
#'   matching \code{state_names}.
#' @param output Optional list defining output conventions. By default,
#'   \code{list(incidence_col = "incidence", cumulative_col = "C")}.
#' @param equations Optional character vector listing model equations for display
#'   in \code{print.epi_model}. When \code{NULL}, the right-hand side function
#'   source is printed instead.
#'
#' @return
#' An object of class \code{"epi_model"}, implemented as a named list.
#'
#' @examples
#' ## Minimal example (structure only)
#' dummy_rhs <- function(time, state, parms) {
#'   list(rep(0, length(state)))
#' }
#'
#' model <- new_epi_model(
#'   name = "DUMMY",
#'   rhs = dummy_rhs,
#'   state_names = c("X"),
#'   par_names = c("alpha")
#' )
#'
#' model
#'
#' @seealso
#' \code{\link{simulate_epi}}
#'
#' @export
new_epi_model <- function(name,
                          rhs,
                          state_names,
                          par_names,
                          lower = NULL,
                          upper = NULL,
                          defaults = NULL,
                          init = NULL,
                          output = list(incidence_col = "incidence"),
                          equations = NULL) {

  stopifnot(is.character(name), length(name) == 1)
  stopifnot(is.function(rhs))
  stopifnot(is.character(state_names), length(state_names) >= 1)
  stopifnot(is.character(par_names), length(par_names) >= 1)

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

  if (!is.null(defaults)) {
    stopifnot(is.numeric(defaults), all(par_names %in% names(defaults)))
    defaults <- defaults[par_names]
  }

  if (!is.null(equations)) {
    stopifnot(is.character(equations))
  }

  structure(
    list(
      name = name,
      rhs = rhs,
      state_names = state_names,
      par_names = par_names,
      lower = lower,
      upper = upper,
      defaults = defaults,
      init = init,
      output = output,
      equations = equations
    ),
    class = "epi_model"
  )
}

#' Print method for epidemic model objects
#' @name print.epi_model
#' @description
#' Provides a concise, human-readable summary of an \code{epi_model} object,
#' including its name, state variables, parameters, and (when available)
#' parameter bounds.
#'
#' This method is automatically called when an \code{epi_model} object is printed
#' at the console.
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
  cat("<epi_model>", x$name, "\n")
  cat("  States: ", paste(x$state_names, collapse = ", "), "\n", sep = "")
  cat("  Params: ", paste(x$par_names, collapse = ", "), "\n", sep = "")

  if (!is.null(x$lower) && !is.null(x$upper)) {
    cat("  Bounds:\n")
    b <- cbind(lower = x$lower, upper = x$upper)
    print(b)
  }

  # --- NUEVO: ecuaciones del modelo (si existen) ---
  if (!is.null(x$equations) && length(x$equations) > 0) {
    cat("  Equations:\n")
    for (ln in x$equations) cat("   ", ln, "\n", sep = "")
  } else {
    # fallback: mostrar el rhs "como código" si no hay ecuaciones guardadas
    cat("  Equations (rhs):\n")
    rhs_txt <- deparse(x$rhs)
    # quita líneas vacías por estética
    rhs_txt <- rhs_txt[nzchar(trimws(rhs_txt))]
    for (ln in rhs_txt) cat("   ", ln, "\n", sep = "")
  }

  # --- Función make_init ---
  if (!is.null(x$make_init) && is.function(x$make_init)) {
    cat("  Initial conditions (make_init):\n")
    mi_txt <- deparse(x$make_init)
    mi_txt <- mi_txt[nzchar(trimws(mi_txt))]
    for (ln in mi_txt) cat("   ", ln, "\n", sep = "")
  } else {
    cat("  Initial conditions: <none>\n")
  }

  invisible(x)
}
