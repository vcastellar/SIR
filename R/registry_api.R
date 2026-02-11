#' Register an epidemiological model
#'
#' Allows users to register custom epi_model objects
#'
#' @param model An object of class "epi_model"
#' @param family Optional family name
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
#' @export
list_models <- function() {
  epi_registry$get_entry_names()
}

#' Get registered model
#' @export
get_model <- function(name) {
  epi_registry$get_entry(name)$model
}
