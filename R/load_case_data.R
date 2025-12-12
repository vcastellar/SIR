#' Cargar datos de casos confirmados
#'
#' Descarga la serie temporal de casos confirmados de COVID-19 y elimina las
#' columnas no necesarias.
#'
#' @param url Dirección del CSV con los datos. Por defecto se usa
#'   `DATA_URL`.
#' @return Un `data.frame` con las columnas de provincia/estado, país y la
#'   serie temporal de casos.
#' @examples
#' datos <- load_case_data()
#' head(datos)
#' @importFrom dplyr select
#' @importFrom utils read.csv
load_case_data <- function(url = DATA_URL) {
  read.csv(url, stringsAsFactors = FALSE) %>%
    select(-Lat, -Long)
}
