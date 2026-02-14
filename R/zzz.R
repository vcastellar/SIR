#' @importFrom stats coef setNames
#' @importFrom graphics plot.new text
#' @importFrom shiny observeEvent
#' @importFrom registry registry

.onLoad <- function(libname, pkgname) {

  register_epi_model(SI_MODEL, family = "SI")
  register_epi_model(SIR_MODEL, family = "SIR")
  register_epi_model(SIRS_MODEL, family = "SIRS")
  register_epi_model(SEIR_MODEL, family = "SEIR")
  register_epi_model(SEIRS_MODEL, family = "SEIRS")
  register_epi_model(SIR_VITAL_MODEL, family = "SIR_VITAL")
}

.register_builtin <- function(model, family) {

  epi_registry$set_entry(
    name = model$name,
    model = model,
    family = family,
    origin = "builtin"
  )
}

