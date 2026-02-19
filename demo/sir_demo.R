###############################################################
# Demo: SIR Model with Vaccination (SIRV)
###############################################################

cat("============================================================\n")
cat("Demo: SIR Model with Vaccination (SIRV)\n")
cat("============================================================\n\n")

cat("This demo illustrates:\n")
cat("- Model construction\n")
cat("- Default parameters\n")
cat("- Incidence derived variable\n")
cat("- Simulation and visualization\n\n")


###############################################################
# Step 1: Model definition
###############################################################

cat("Step 1: Defining the model structure\n")
cat("------------------------------------------------------------\n")
cat("Compartments:\n")
cat("S = Susceptible\n")
cat("I = Infectious\n")
cat("R = Recovered\n")
cat("V = Vaccinated\n\n")

cat("Derived variable:\n")
cat("incidence = beta * S * I / N\n")
cat("This represents the number of new infections per unit time.\n\n")


sirv_rhs <- function(time, state, parms) {

  with(as.list(c(state, parms)), {

    N <- S + I + R + V

    incidence <- beta * S * I / N

    dS <- -incidence - nu * S
    dI <-  incidence - gamma * I
    dR <-  gamma * I
    dV <-  nu * S

    list(
      c(dS, dI, dR, dV),
      incidence = incidence
    )
  })
}

SIRV_MODEL <- epi_model(
  name      = "SIRV",
  rhs       = sirv_rhs,
  par_names = c("beta", "gamma", "nu"),
  states    = c("S", "I", "R", "V"),
  derived   = "incidence",
  defaults  = c(beta = 0.3, gamma = 0.1, nu = 0.01)
)

cat("Model successfully created.\n\n")


###############################################################
# Step 2: Default parameters
###############################################################

cat("Step 2: Understanding default parameters\n")
cat("------------------------------------------------------------\n")
cat("Defaults are baseline parameter values stored inside\n")
cat("the model object.\n")
cat("They represent a reasonable epidemiological scenario\n")
cat("and allow simulation without manually specifying parameters.\n\n")

print(SIRV_MODEL$defaults)
cat("\n")


###############################################################
# Step 3: Simulation
###############################################################

cat("Step 3: Running a simulation\n")
cat("------------------------------------------------------------\n")

init <- c(S = 1e5, I = 10, R = 0, V = 0)

cat("Initial conditions:\n")
print(init)
cat("\n")

sim <- simulate_epi(
  model = SIRV_MODEL,
  times = 0:200,
  parms = SIRV_MODEL$defaults,
  init  = init
)

cat("Simulation completed.\n\n")


###############################################################
# Step 4: Outputs
###############################################################

cat("Step 4: Inspecting outputs\n")
cat("------------------------------------------------------------\n")
cat("The simulation includes:\n")
cat("- State variables (S, I, R, V)\n")
cat("- Derived variable: incidence\n\n")

print(summary(sim))
cat("\n")


###############################################################
# Step 5: Plot
###############################################################

cat("Step 5: Plotting trajectories\n")
cat("------------------------------------------------------------\n\n")

plot(sim)

cat("\nThe incidence curve corresponds to new infections.\n")
cat("It represents the transition from S to I.\n\n")



###############################################################
# Register model in the internal registry
###############################################################

register_epi_model(SIRV_MODEL)
list_models()
cat("Model registered in internal registry.\n\n")

run_epi_app()


