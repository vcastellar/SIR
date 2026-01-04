
suppressPackageStartupMessages({
  library(deSolve)
})

# helper: operador "si NULL entonces"
`%||%` <- function(a, b) if (is.null(a)) b else a

# library(dplyr)
# infectedCV <- read.csv("data/infectedCV.csv")
# x <- infectedCV %>% group_by(fecha) %>% summarise(N = sum(CasosAcum))
# x <- diff(x$N)[1:350]
# plot(x)
#x <- simulate_sir()$incidence_obs$inc
# saveRDS(x, file = "./data/incidencias.rds")
