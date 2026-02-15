#' Register an epidemiological model
#'
#' @param model An object of class "epi_model".
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



#' List registered models
#'
#' @return Character vector with registered model names.
#' @export
list_models <- function() {
  ls(envir = .epi_registry)
}


#' Get a registered model
#'
#' @param name Character scalar with model name.
#' @return An object of class "epi_model".
#' @export
get_model <- function(name) {

  if (!exists(name, envir = .epi_registry, inherits = FALSE)) {
    stop("Model not found in registry.")
  }

  get(name, envir = .epi_registry)
}


