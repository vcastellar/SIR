
suppressPackageStartupMessages({
  library(deSolve)
})

# helper: operador "si NULL entonces"
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Safe maximum ignoring non-finite values
#'
#' Internal helper to compute a maximum without producing warnings when
#' the input vector is empty or contains only non-finite values.
#'
#' @keywords internal
#' @noRd
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}



#' Safe minimum ignoring non-finite values
#'
#' Internal helper to compute a minimum without producing warnings when
#' the input vector is empty or contains only non-finite values.
#'
#' @keywords internal
#' @noRd
safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

# library(dplyr)
# infectedCV <- read.csv("data/infectedCV.csv")
# x <- infectedCV %>% group_by(fecha) %>% summarise(N = sum(CasosAcum))
# x <- diff(x$N)[1:350]
# plot(x)
#x <- simulate_sir()$incidence_obs$inc
# saveRDS(x, file = "./data/incidencias.rds")
