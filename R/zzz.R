#' @importFrom stats coef setNames
#' @importFrom graphics plot.new text
#' @importFrom shiny observeEvent

.onLoad <- function(libname, pkgname) {

  register_epi_model(SI_MODEL)
  register_epi_model(SIR_MODEL)
  register_epi_model(SIRS_MODEL)
  register_epi_model(SEIR_MODEL)
  register_epi_model(SEIRS_MODEL)
  register_epi_model(SIR_VITAL_MODEL)
}
