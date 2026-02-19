#' Return built-in epidemiological models
#'
#' @keywords internal
#' @export
.get_builtin_models <- function() {

  list(
    SI    = SI_MODEL,
    SIR   = SIR_MODEL,
    SIRS  = SIRS_MODEL,
    SEIR  = SEIR_MODEL,
    SEIRS = SEIRS_MODEL,
    SIR_VITAL = SIR_VITAL_MODEL
  )
}



#' example of new epi model
#' @keywords internal
#' @examples
#' ## ---------------------------------------------------------------
#' ## Define and register a new epidemiological model (SIRD)
#' ## ---------------------------------------------------------------
#'
#' sird_rhs <- function(time, state, parms) {
#'   with(as.list(c(state, parms)), {
#'     N <- S + I + R
#'     lambda <- beta * S * I / N
#'
#'     dS <- -lambda
#'     dI <-  lambda - gamma * I - mu * I
#'     dR <-  gamma * I
#'     dD <-  mu * I
#'
#'     list(
#'       c(dS, dI, dR, dD),
#'       incidence = lambda
#'     )
#'   })
#' }
#'
#' SIRD_MODEL <- epi_model(
#'   name      = "SIRD",
#'   rhs       = sird_rhs,
#'   states    = c("S", "I", "R", "D"),
#'   derived   = c("incidence"),
#'   par_names = c("beta", "gamma", "mu"),
#'   defaults  = c(beta = 0.3, gamma = 0.1, mu = 0.01),
#'   init      = c(S = 1e6, I = 10, R = 0, D = 0)
#' )
#'
#' ## Simulate using explicit initial conditions
#' sim <- simulate_epi(
#'   model = SIRD_MODEL,
#'   times = 0:200,
#'   parms = c(beta = 0.25, gamma = 0.1, mu = 0.02),
#'   init  = SIRD_MODEL$init
#' )
#'
#' plot(sim)
#' plot(sim, what = "derived")
#' run_epi_app(list(SIRD = SIRD_MODEL))



