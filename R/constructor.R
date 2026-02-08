#-------------------------------------------------------------------------------
# Epidemiological role vocabulary (state roles only)
#-------------------------------------------------------------------------------
#' Epidemiological role vocabulary
#'
#' Defines the internal vocabulary of epidemiological roles supported by the
#' framework.
#'
#' This vocabulary is restricted to *state roles only*. Epidemiological roles
#' describe the semantic meaning of compartmental state variables (e.g.
#' susceptible, infectious) and are used to provide a model-independent way of
#' referring to core epidemiological concepts.
#'
#' Flows are intentionally excluded from the role system. Flow variables are
#' treated as implementation details or optional model outputs and are not
#' required to have epidemiological roles.
#'
#' @details
#' ## State roles
#'
#' The vocabulary defines the set of supported epidemiological state roles and
#' their default associated state variable names. The returned vector maps role
#' names (keys) to default state names (values).
#'
#' The following state roles are currently supported:
#'
#' \describe{
#'   \item{\code{susceptible}}{
#'     Individuals who are not infected and are at risk of infection.
#'     By default, associated with the state variable \code{"S"}.
#'   }
#'   \item{\code{exposed}}{
#'     Individuals who have been infected but are not yet infectious.
#'     By default, associated with the state variable \code{"E"}.
#'   }
#'   \item{\code{infectious}}{
#'     Individuals who are currently infectious and capable of transmitting
#'     the pathogen.
#'     By default, associated with the state variable \code{"I"}.
#'   }
#'   \item{\code{recovered}}{
#'     Individuals who have recovered from infection and are no longer
#'     infectious.
#'     By default, associated with the state variable \code{"R"}.
#'   }
#'   \item{\code{deceased}}{
#'     Individuals who have died as a consequence of the disease.
#'     By default, associated with the state variable \code{"D"}.
#'   }
#' }
#'
#' ## Use in \code{epi_model()}
#'
#' The vocabulary is used by \code{epi_model()} to:
#'
#' \itemize{
#'   \item Infer state roles automatically when \code{roles = NULL}.
#'   \item Validate user-defined state role assignments.
#'   \item Ensure that all model states have an epidemiological role.
#' }
#'
#' All state variables declared in an \code{epi_model} must be identifiable
#' through this vocabulary or explicitly assigned by the user. If a state
#' cannot be mapped to a role, model construction fails.
#'
#' @return
#' A named character vector mapping epidemiological state roles (names) to
#' default state variable names (values).
#'
#' @seealso
#' \code{\link{epi_model}}
#'
#' @keywords internal


epi_role_vocab <- function() {
  c(
    susceptible = "S",
    exposed     = "E",
    infectious  = "I",
    recovered   = "R",
    deceased    = "D"
  )
}


#-------------------------------------------------------------------------------
# Epidemiological model constructor
#-------------------------------------------------------------------------------
#' Construct an epidemiological model
#'
#' Creates an \code{epi_model} object describing the structure of a deterministic
#' compartmental epidemiological model. The model definition includes the
#' system of equations, state variables, optional flow variables, parameter
#' names, and a semantic mapping between state variables and epidemiological
#' roles.
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
#'   used internally by the model or returned as outputs. Flow variables are
#'   not required to have epidemiological roles.
#'
#' @param roles Named list mapping *epidemiological roles* to state variables.
#'   Names correspond to epidemiological state roles (see
#'   \code{\link{epi_role_vocab}}), and values correspond to names of variables
#'   defined in \code{states}.
#'
#'   If \code{roles = NULL} (default), roles are inferred automatically using
#'   the internal epidemiological role vocabulary. In all cases, every state
#'   variable must be associated with exactly one epidemiological role.
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
#' @param outputs Optional character vector giving the names of model variables
#'   (states and/or flows) to be returned by default in simulations. If
#'   \code{NULL}, all states and flows are included.
#'
#' @details
#' ## Epidemiological roles
#'
#' Epidemiological roles provide a semantic layer on top of the model states,
#' allowing downstream tools (metrics, summaries, plotting functions) to refer
#' to core epidemiological concepts (e.g. susceptible, infectious) in a
#' model-independent way.
#'
#' Roles apply \strong{exclusively to state variables}. Flow variables are not
#' part of the role system and are treated as implementation details or optional
#' model outputs.
#'
#' The set of supported epidemiological roles and their default associated state
#' names are defined by \code{\link{epi_role_vocab}}.
#'
#' ### Automatic role assignment
#'
#' When \code{roles = NULL}, the constructor attempts to assign roles
#' automatically by matching the default state names defined in
#' \code{epi_role_vocab()} against the declared \code{states}.
#'
#' Model construction fails if any state variable cannot be associated with an
#' epidemiological role.
#'
#' ### User-defined roles
#'
#' When \code{roles} is provided by the user:
#'
#' \itemize{
#'   \item User-defined role assignments take precedence over automatic
#'   inference.
#'   \item Role names must belong to the epidemiological role vocabulary.
#'   \item Role targets must correspond to state variables declared in
#'   \code{states}.
#' }
#'
#' Any roles not explicitly defined by the user are completed automatically
#' using the vocabulary. All states must be identifiable after completion.
#'
#' @return
#' An object of class \code{"epi_model"} containing the model definition.
#'
#' @seealso
#' \code{\link{simulate_epi}},
#' \code{\link{epi_role_vocab}}
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
                      init = NULL,
                      outputs = NULL) {

  ## --- basic checks ----------------------------------------------------------
  stopifnot(is.character(name), length(name) == 1)
  stopifnot(is.function(rhs))
  stopifnot(is.character(par_names), length(par_names) >= 1)
  stopifnot(is.character(states), length(states) >= 1)
  stopifnot(length(unique(states)) == length(states))

  if (missing(flows) || length(flows) == 0) {
    flows <- character(0)
  }

  ## a variable cannot be both state and flow
  stopifnot(
    is.character(flows),
    length(unique(flows)) == length(flows),
    !any(flows %in% states)
  )

  #-----------------------------------------------------------------------------
  ## Roles handling (STATES ONLY)
  #-----------------------------------------------------------------------------

  vocab <- epi_role_vocab()
  model_states <- states

  ## -------------------------------------------------------------------------
  ## Case 1: roles = NULL  → infer from vocabulary
  ## -------------------------------------------------------------------------
  if (is.null(roles)) {

    roles <- list()

    for (role in names(vocab)) {
      state <- vocab[[role]]

      if (state %in% model_states) {
        roles[[role]] <- state
      }
    }

    ## all states must have a role
    missing_states <- setdiff(model_states, unlist(roles))

    if (length(missing_states) > 0) {
      stop(
        "Unable to assign epidemiological roles to the following states: ",
        paste(missing_states, collapse = ", ")
      )
    }

    ## -------------------------------------------------------------------------
    ## Case 2 & 3: roles defined by the user (partial or complete)
    ## -------------------------------------------------------------------------
  } else {

    ## --- structure ----------------------------------------------------------
    if (!is.list(roles) || is.null(names(roles))) {
      stop("`roles` must be a named list.")
    }

    ## --- role names must exist in vocabulary --------------------------------
    unknown_roles <- setdiff(names(roles), names(vocab))
    if (length(unknown_roles) > 0) {
      stop(
        "Unknown epidemiological roles: ",
        paste(unknown_roles, collapse = ", ")
      )
    }

    ## --- role targets must be states ----------------------------------------
    unknown_states <- setdiff(unlist(roles), model_states)
    if (length(unknown_states) > 0) {
      stop(
        "Roles refer to unknown states: ",
        paste(unknown_states, collapse = ", ")
      )
    }

    ## --- complete missing roles using vocabulary ----------------------------
    missing_roles <- setdiff(names(vocab), names(roles))

    for (role in missing_roles) {
      state <- vocab[[role]]

      if (state %in% model_states) {
        roles[[role]] <- state
      }
    }

    ## --- all states must have a role ----------------------------------------
    missing_states <- setdiff(model_states, unlist(roles))

    if (length(missing_states) > 0) {
      stop(
        "Unable to assign epidemiological roles to the following states: ",
        paste(missing_states, collapse = ", ")
      )
    }
  }

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

  ## --- outputs ---------------------------------------------------------------
  if (is.null(outputs)) {
    outputs <- unique(c(states, flows))
  }
  stopifnot(is.character(outputs), all(states %in% outputs))

  structure(
    list(
      name = name,
      rhs = rhs,
      par_names = par_names,
      states = states,
      flows = flows,
      outputs = outputs,
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
