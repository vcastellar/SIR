
suppressPackageStartupMessages({
  library(deSolve)
})

infected <- readRDS("./data/casos_ts.rds")
infected <- infected[100:360]
plot(diff(infected), type = 'l')
infected <- round(infected)
# saveRDS(infected, file = "./data/casos_ts.rds")
