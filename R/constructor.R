## ============================================================================
## Epidemiological roles vocabulary
## ============================================================================

#' Epidemiological roles vocabulary
#'
#' @description
#' Defines the standard vocabulary of epidemiological roles supported by the
#' package. Roles assign semantic epidemiological meaning to **state variables**
#' and **flows** of an epidemic model (e.g. susceptible, infectious, incidence).
#'
#' Roles are used to ensure semantic consistency across models and to enable
#' generic epidemiological metrics and visualisations.
#'
#' @details
#' Roles can only be assigned to:
#' \itemize{
#'   \item state variables (e.g. \code{S}, \code{I}, \code{R}),
#'   \item flows or auxiliary outputs (e.g. incidence, deaths).
#' }
#'
#' Epidemiological roles **must not** be assigned to model parameters
#' (e.g. transmission or recovery rates). Models declaring roles outside this
#' vocabulary, or assigning roles to parameters, will be rejected by
#' \code{epi_model()}.
#'
#' The following roles are currently supported:
#'
#' \itemize{
#'   \item \strong{susceptible}: Individuals susceptible to infection.
#'   \item \strong{exposed}: Infected but not yet infectious individuals.
#'   \item \strong{infectious}: Infectious individuals.
#'   \item \strong{recovered}: Individuals recovered with immunity (temporary or permanent).
#'   \item \strong{deceased}: Individuals who died due to the disease.
#'   \item \strong{incidence}: Flow representing new infections.
#'   \item \strong{deaths}: Flow representing disease-induced deaths.
#' }
#'
#' This vocabulary may be extended in future versions of the package.
#'
#' @return
#' A character vector of valid epidemiological role names.
#'
#' @keywords internal

.epi_role_vocab <- function() {
  c(
    "susceptible",
    "exposed",
    "infectious",
    "recovered",
    "deceased",
    "incidence",
    "deaths"
  )
}


## ============================================================================
## Epidemic model constructor
## ============================================================================

#' Define an epidemic model
#'
#' @description
#' Defines a deterministic compartmental epidemic model governed by a system
#' of ordinary differential equations (ODEs). The model explicitly declares its
#' state variables, parameters, model outputs, and optional epidemiological roles.
#'
#' Epidemiological roles assign semantic meaning (e.g. infectious, incidence)
#' to specific state variables or auxiliary outputs, enabling generic
#' epidemiological metrics and summaries to operate consistently across models.
#'
#' @param name Character scalar. Human-readable name of the model (e.g. "SIR").
#' @param rhs Function defining the ODE system. Must have signature
#'   \code{function(time, state, parms)} and return a list whose first element
#'   is the vector of state derivatives.
#' @param states Character vector of state variable names.
#' @param par_names Character vector of parameter names.
#' @param outputs Character vector of named outputs returned by the RHS.
#'   Must include all state variables.
#' @param roles Optional named list mapping epidemiological roles to variable
#'   names (states or outputs). Role names must belong to the standard
#'   vocabulary defined in \code{.epi_role_vocab()}.
#' @param lower Optional named numeric vector of lower parameter bounds.
#' @param upper Optional named numeric vector of upper parameter bounds.
#' @param defaults Optional named numeric vector of default parameter values.
#' @param init Optional named numeric vector of default initial conditions.
#'
#' @details
#' ## Epidemiological roles
#' Roles provide a semantic layer on top of model variables. For example,
#' different models may use different variable names for infectious individuals,
#' but assigning the \code{"infectious"} role allows epidemiological metrics to
#' remain model-agnostic.
#'
#' Supported roles are:
#' \itemize{
#'   \item \strong{susceptible}
#'   \item \strong{exposed}
#'   \item \strong{infectious}
#'   \item \strong{recovered}
#'   \item \strong{deceased}
#'   \item \strong{incidence}
#'   \item \strong{deaths}
#' }
#'
#' Not all roles must be defined for every model. Metrics depending on missing
#' roles should fail explicitly.
#'
#' @return
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' sir_model <- epi_model(
#'   name = "SIR",
#'   rhs  = sir_rhs,
#'   states = c("S", "I", "R"),
#'   par_names   = c("beta", "gamma"),
#'   outputs     = c("S", "I", "R", "incidence"),
#'   roles = list(
#'     susceptible = "S",
#'     infectious  = "I",
#'     recovered   = "R",
#'     incidence   = "incidence"
#'   )
#' )
#'
#' @export
epi_model <- function(name,
                      rhs,
                      par_names,
                      states,
                      flows = character(0),
                      roles = NULL,
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

  if (missing(flows) || length(flows) == 0) {
    flows <- character(0)
  }

  # una variable no puede ser state y flow a la vez
  stopifnot(
    is.character(flows),
    length(unique(flows)) == length(flows),
    !any(flows %in% states)
  )

  ## --- roles validation ------------------------------------------------------
  if (!is.null(roles)) {

    if (!is.list(roles) || is.null(names(roles))) {
      stop("`roles` must be a named list.")
    }

    valid_roles <- .epi_role_vocab()
    role_names  <- names(roles)

    unknown_roles <- setdiff(role_names, valid_roles)
    if (length(unknown_roles) > 0) {
      stop(
        "Unknown epidemiological roles: ",
        paste(unknown_roles, collapse = ", "),
        ". Valid roles are: ",
        paste(valid_roles, collapse = ", ")
      )
    }

    role_vars <- unlist(roles, use.names = FALSE)

    unknown_vars <- setdiff(role_vars, c(states, flows))
    if (length(unknown_vars) > 0) {
      stop(
        "Roles refer to unknown outputs: ",
        paste(unknown_vars, collapse = ", ")
      )
    }

    if (any(duplicated(role_vars))) {
      stop("Each role must map to a unique variable.")
    }
  }

  ## --- parameter bounds ------------------------------------------------------
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
      flows = flows,
      roles = roles,
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
#' the declared model outputs, and the underlying system of differential
#' equations.
#'
#' This method is automatically called when an object of class
#' \code{"epi_model"} is printed at the console.
#'
#' @details
#' The \code{epi_model} class explicitly declares the set of model outputs via
#' the \code{outputs} field. These outputs may include state variables, derived
#' quantities such as incidence, or any other observable returned by the
#' model's right-hand side (\code{rhs}) function.
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

  ## --- outputs ---------------------------------------------------------------
  if (!is.null(x$outputs)) {
    cat("  Outputs:  ", paste(x$outputs, collapse = ", "), "\n", sep = "")
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
