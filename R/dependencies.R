#' Dependencias del paquete
#'
#' Carga silenciosamente las dependencias necesarias para ejecutar el modelo
#' SIR y la curva logística.
#' @keywords internal
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(reshape2)
  library(xts)
  library(zoo)
  library(deSolve)
})

#' URL de los datos públicos
#'
#' Fuente base con los casos confirmados de COVID-19.
#' @keywords internal
DATA_URL <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
