#' Register an epidemiological model
#'
#' @description
#' Allows users to register custom \'epi_model\' objects in the package
#' registry so they can be retrieved and used by other utilities.
#'
#' @param model An object of class \code{"epi_model"}.
#' @param family Optional character scalar indicating a model family.
#'
#' @return
#' Invisibly returns \code{TRUE} when registration succeeds.
#' @export
register_epi_model <- function(model, family = NULL) {

  stopifnot(inherits(model, "epi_model"))

  name <- model$name

  if (epi_registry$has_entry(name)) {
    stop("A model with this name is already registered.")
  }

  epi_registry$set_entry(
    name = name,
    model = model,
    family = family
  )

  invisible(TRUE)
}

#' Remove a registered epidemiological model
#'
#' @description
#' Removes a model from the registry by name. Built-in models cannot
#' be removed.
#'
#' @param name Character scalar with the model name to remove.
#'
#' @return
#' Invisibly returns \code{TRUE} when removal succeeds.
#' @export
remove_epi_model <- function(name) {

  if (!epi_registry$has_entry(name)) {
    stop("Model not found.")
  }

  entry <- epi_registry$get_entry(name)

  if (entry$origin == "builtin") {
    stop("Built-in models cannot be removed.")
  }

  epi_registry$delete_entry(name)

  invisible(TRUE)
}



#' List registered models
#'
#' @description
#' Returns the names of all models currently stored in the registry.
#'
#' @return
#' Character vector with registered model names.
#' @export
list_models <- function() {
  epi_registry$get_entry_names()
}

#' Get registered model
#'
#' @description
#' Retrieves a model from the registry by name.
#'
#' @param name Character scalar with the model name to retrieve.
#'
#' @return
#' An object of class \code{"epi_model"}.
#' @export
get_model <- function(name) {
  epi_registry$get_entry(name)$model
}
