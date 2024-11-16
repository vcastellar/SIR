
library(dplyr)
library(lubridate)
library(tidyr)
library(reshape2)
library(xts)
deaths <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv",
                   stringsAsFactors = FALSE)
head(deaths)
deaths$Lat <- NULL
deaths$Long <- NULL


deathsLong <- melt(deaths, id.vars = c("Province.State", "Country.Region"))
deathsLong$variable <- gsub("^X", "", deathsLong$variable)
deathsLong$variable <- gsub("\\.", "-", deathsLong$variable)
deathsLong$variable <- gsub("-20$", "-2020", deathsLong$variable)

deathsLong$variable <- as.Date(as.character(deathsLong$variable), format = "%m-%d-%Y")

agregados  <- deathsLong %>% 
  
  group_by(Country.Region, variable) %>% 
  summarise(deaths = sum(value)) %>% 
  as.data.frame

sort(table(agregados$Country.Region))

n_init_infect <- 100

kk <- agregados %>% group_by(Country.Region) %>% 
  filter(deaths > n_init_infect) %>% 
  summarise(alpha = sum(diff(log(tail(deaths, 4))), na.rm = TRUE)/3) %>% data.frame() %>% 
  mutate(z = log(2) / log(1 + alpha)) %>% 
  arrange(z)


espanya <- deathsLong %>% filter(Country.Region == "Spain")
N <- 47e6

Infected <- round(espanya$value)
Infected <- Infected[Infected > n_init_infect]


Fecha <- as.Date(espanya$variable[which(espanya$value > n_init_infect)[1]])
Fechas_Obs <- seq.Date(from = Fecha, length.out = length(Infected), by = "day")

serie <- ts(Infected, start = c(2020, as.numeric(format(Fechas_Obs[1], "%j"))),
            frequency = 365)
plot(serie)

serie <- xts(x = Infected, order.by = Fechas_Obs)
plot(serie)
plot(log2(serie))
as.numeric(diff(log(serie)))
plot(diff(log(serie)))

IncPorc <- rollapply(Infected, width = 2, FUN = function(x) (x[2] - x[1])/x[1], align = "left")
alpha <- rollapply(IncPorc, width = 4, FUN = mean)
DiasDobla <- log(2) / (log(1 + alpha))
serieDobla <- xts(x = DiasDobla, order.by = Fechas_Obs[-(1:4)])
serieDobla[is.infinite(serieDobla)] <- NA
plot(serieDobla, main = "nº de días estimados en lo que se dobla el nº de infectados", ylim = c(0, 1e3))
# diás estimados en los que se doblan el número de infectados
log(2) / log(1 + IncPorc)
barplot(diff(serie), las = 2, cex.names = 0.7)

# Infected <- c(45, 62, 121, 198, 291, 440, 571, 830, 1287, 1975, 2744, 4515, 5974, 7711, 9692, 11791, 14380, 17205, 20440)

# - S: número de susceptibles
# - I: número de infectados
# - R: número de recobrados
Day <- 1:(length(Infected))
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}




library(deSolve)
init <- c(S = N-Infected[1], I = Infected[1], R = 0)
RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, 
             times = Day, 
             func = SIR, 
             parms = parameters)
  fit <- out[ , 3]
  mean((Infected - fit)^2)
}

# optimize with some sensible conditions
Opt <- optim(c(0.5, 0.5), 
             RSS, 
             method = "Nelder-Mead",
             control = list(maxit = 500, factr = 1e4)) 
Opt$message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par
##      beta     gamma 
## 0.6746089 0.3253912
## R0
Opt_par[1] / Opt_par[2]

t <- 1:(90 + length(Infected)) # time in days
fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))



Fechas <- seq.Date(from = Fecha, length.out = nrow(fit), by = "day")
fit$time <- Fechas


plot(fit$time, fit$S, type = "l", lwd = 2, ylim = c(0, N))
lines(fit$time, fit$I, type = "l", lwd = 2, col = "red")
lines(fit$time, fit$R, type = "l", lwd = 2, col = "green")
points(Fechas[1:length(Infected)], Infected)
grid()
legend("left", c("Susceptibles", "Infectados", "Recobrados"), lty = 1, lwd = 2, col = 1:3, inset = 0.05)

plot(fit$time, fit$I, type = "l", lwd = 2, col = "red", pch = 19)
points(Fechas[1:length(Infected)], Infected)
grid()



logistic <- function(parms, x) {
  parms[3]/(1 + exp(- (x - parms[2]) / parms[1]))
}

parms <- c(A = 1, B = 1, C = 2e6)
RSS <- function(parms) {
  out <- logistic(parms, x = Day)
  fit <- out
  sum((Infected - fit)^2)
}

Opt <- optim(c(1, 1, 2e6), RSS, method = "L-BFGS-B") 

fit <- logistic(Opt$par, 1:length(Fechas))
plot(as.Date(Fechas), fit)
points(Fechas_Obs, Infected, pch = 19, col = "red")
grid()







Fechas <- seq.Date(from = Fecha, length.out = length(Infected), by = "day")
points(Fechas, Infected)
legend("right", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col, inset = 0.05)
title("SIR model 2019-nCoV Espa?a", outer = TRUE, line = -2)
