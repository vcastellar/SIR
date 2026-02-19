#' Register an epidemiological model
#'
#' @description
#' Registers an object of class \code{"epi_model"} in the internal
#' model registry of the package.
#'
#' The registry is an in-memory environment that stores epidemiological
#' models by name. Once registered, a model becomes discoverable through
#' \code{\link{list_models}} and retrievable via \code{\link{get_model}}.
#'
#' This mechanism allows users to extend the framework with custom models
#' that can be used in simulation workflows or in the Shiny application.
#'
#' @param model An object of class \code{"epi_model"}.
#'
#' @details
#' The model is stored using its \code{name} field as identifier.
#' Model names must be unique within the current R session.
#'
#' Attempting to register a model with a name that already exists in the
#' registry results in an error.
#'
#' The registry is session-scoped and is reinitialized when the package
#' is reloaded.
#'
#' @return
#' Invisibly returns \code{TRUE} on success.
#'
#' @examples
#' \dontrun{
#' m <- SIR_MODEL
#' m$name <- "SIR_custom"
#'
#' register_epi_model(m)
#' list_models()
#' }
#'
#' @seealso
#' \code{\link{list_models}},
#' \code{\link{get_model}},
#' \code{\link{epi_model}}
#'
#' @export
register_epi_model <- function(model) {

  stopifnot(inherits(model, "epi_model"))

  name <- model$name

  if (exists(name, envir = .epi_registry, inherits = FALSE)) {
    stop("A model with this name is already registered.")
  }

  assign(name, model, envir = .epi_registry)

  invisible(TRUE)
}



#' List registered epidemiological models
#'
#' @description
#' Returns the names of all epidemiological models currently stored
#' in the internal registry.
#'
#' @details
#' The registry contains:
#' \itemize{
#'   \item Built-in models registered at package load time.
#'   \item User-defined models registered via \code{\link{register_epi_model}}.
#' }
#'
#' The returned value is a character vector containing model names.
#'
#' @return
#' A character vector with the names of registered models.
#'
#' @examples
#' \dontrun{
#' library(SIR)
#' list_models()
#' }
#'
#' @seealso
#' \code{\link{register_epi_model}},
#' \code{\link{get_model}}
#'
#' @export
list_models <- function() {
  ls(envir = .epi_registry)
}


#' Retrieve a registered epidemiological model
#'
#' @description
#' Retrieves a previously registered epidemiological model from the
#' internal registry by name.
#'
#' @param name Character scalar specifying the model name.
#'
#' @details
#' The model must have been registered either at package load time
#' (built-in models) or via \code{\link{register_epi_model}}.
#'
#' If no model with the specified name exists in the registry,
#' an error is raised.
#'
#' @return
#' An object of class \code{"epi_model"}.
#'
#' @examples
#' \dontrun{
#' model <- get_model("SIR")
#' print(model)
#' }
#'
#' @seealso
#' \code{\link{register_epi_model}},
#' \code{\link{list_models}}
#'
#' @export
get_model <- function(name) {

  if (!exists(name, envir = .epi_registry, inherits = FALSE)) {
    stop("Model not found in registry.")
  }

  get(name, envir = .epi_registry)
}


#' Remove a registered epidemiological model
#'
#' @description
#' Removes a previously registered epidemiological model from the
#' internal registry.
#'
#' @param name Character scalar specifying the model name to remove.
#'
#' @details
#' The function deletes the model from the internal in-memory registry.
#' This affects only the current R session.
#'
#' Attempting to remove a model that does not exist results in an error.
#'
#' Built-in models can also be removed, but they will be restored
#' the next time the package is reloaded.
#'
#' @return
#' Invisibly returns \code{TRUE} on success.
#'
#' @examples
#' \dontrun{
#' register_epi_model(SIR_MODEL)
#' list_models()
#'
#' unregister_epi_model("SIR")
#' list_models()
#' }
#'
#' @seealso
#' \code{\link{register_epi_model}},
#' \code{\link{list_models}},
#' \code{\link{get_model}}
#'
#' @export
unregister_epi_model <- function(name) {

  if (!exists(name, envir = .epi_registry, inherits = FALSE)) {
    stop("Model not found in registry.")
  }

  rm(list = name, envir = .epi_registry)

  invisible(TRUE)
}




