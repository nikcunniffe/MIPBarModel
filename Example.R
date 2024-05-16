####
#### Script to run MIP, BAR and MIP-BAR models under default parameterisation
#### 

#### This script reproduces Figure 4a. It requires packages ggplot2 and deSolve,
#### which can be installed with the command: install.packages(c("ggplot2", "deSolve"))


rm(list=ls())

# packages
library(deSolve)
library(ggplot2)

#### Define functions for RHS of model ODES ####

MIP_ode <- function(times, states, par) {
  
  ### Define terms for dI/dt
  S <- par[["H"]] - states[["I"]] # number of susceptible plants
  
  plant_become_infected <- par[["phi"]] * par[["b"]] * states[["Z"]] * S /
    (S + par[["v"]]*states[["I"]])
  
  plant_death <- par[["gamma"]] * states[["I"]]

  # dI/dt (Plants state equation)
  dI <- plant_become_infected - plant_death
  
  ### Define terms for dZ/dt
  X <- par[["A"]] - states[["Z"]] # number of un-infective aphids
  
  aphid_virus_acquisition <- par[["phi"]] * par[["a"]] * 
    (1 - par[["epsilon"]]*par[["omega"]]) * X * par[["v"]]*states[["I"]] /
    (S + par[["v"]]*states[["I"]])
  
  aphid_virus_loss <- par[["tau"]] * states[["Z"]]
  
  # dZ/dt (Aphids state equation)
  dZ <- aphid_virus_acquisition - aphid_virus_loss
  
  return(list(c(dI, dZ)))
  
}

BAR_ode <- function(times, states, par) {
  
  # STATES
  I <- states[["I"]]
  
  # PARAMETERS
  gamma <- par[["gamma"]]
  theta <- par[["theta"]]
  H <- par[["H"]]
  v <- par[["v"]]
  epsilon <- par[["epsilon"]]
  omega <- par[["omega"]]
  a <- par[["a"]]
  b <- par[["b"]]
  A <- par[["A"]]
  
  # DEFINE PARAMETERS NEEDED FOR STATE EQUATIONS
  # I_hat = weighted number of infected plants, accounting for attraction towards infected plants
  i_hat <- v*I / (H - I + v*I)
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly et al. 2019, Appendix S1)
  xI <- (a*b*i_hat*(1 - epsilon*omega)*(1 - i_hat))/(omega*(1 - i_hat*(1 - epsilon)))
  
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A)*xI - gamma*I
  
  return(list(di))
}

MIP_BAR_ode <- function(times, states, par) {
  
  omega <- par[["omega"]]
  eta <- par[["eta"]]
  H <- par[["H"]]
  A <- par[["A"]]
  b <- par[["b"]]
  a <- par[["a"]]
  gamma <- par[["gamma"]]
  v <- par[["v"]]
  epsilon <- par[["epsilon"]]
  rho <- par[["rho"]] # probability of losing infectivity after probing
  
  # states
  I <- states[["I"]]
  Z <- states[["Z"]]
  
  S <- H - I # number of susceptible plants
  X <- A - Z # number of un-infective aphids
  
  # expression for rate of aphid diserpsals
  phi_hat <- (S + v*I) / (omega*eta*(S + v*epsilon*I))
  
  # expression for rate of aphid infectivity loss
  tau_hat <- phi_hat*((rho*(1 - a*(1 - epsilon*omega)) + (1 - rho)*epsilon*omega)*v*I +
                        (rho + (1 - rho)*omega)*S)/(S + v*I)
  
  dI <- phi_hat*Z*b*S/(S + v*I) - gamma*I
  dZ <- phi_hat*X*a*(1 - epsilon*omega)*v*I/(S + v*I) - tau_hat*Z
    
  
  return(list(c(dI, dZ)))
}


#### Define parameters, initial states and timeframe ####

## Parameters
parms <- c(
  
  # POPULATION SIZES
  H = 10000, # plants
  A = 500, # aphids
  
  # VIRUS ACQUISITION/INOCULATION PROBABILITY
  b = 0.5, 
  a = 0.5,
  
  # PLANT DEATH/REPLANTING RATE
  gamma = 0.03,
  
  # VECTOR PREFERENCE (ATTRACTIVENESS + ACCEPTABILITY)
  v = 1,
  epsilon = 1,
  
  # APHID DISPERSAL
  eta = 0.8333, # aphid feed length (MIP-BAR model)
  phi = 1/(0.2*0.8333), # aphid dispersal rate (MIP model)
  theta = 1/0.8333, # aphid initiation of feeding dispersal rate (BAR model)
  
  # APHID INFECTIVITY LOSS
  tau = 4.8, # aphid infectivity loss rate (MIP model)
  rho = 1, # probability of infectivity loss from probing (MIP-BAR model)
  omega = 0.2 # probability of infectivity loss from feeding
)

## Initial states
init_I <- 1
init_Z <- 0

init_states <- c(I = init_I, Z = init_Z)
init_states_BAR <- c(I = init_I)

## Timeframe
tf <- 600
times <- seq(0, tf, length.out = 200)

#### Run models ####

MIP_trajec <- data.frame(deSolve::ode(y = init_states,
                                      times = times,
                                      func = MIP_ode,
                                      parms = parms))
BAR_trajec <- data.frame(deSolve::ode(y = init_states_BAR,
                                      times = times,
                                      func = BAR_ode,
                                      parms = parms))
BAR_trajec$Z <- NA # add Z column filled with NAs to make combining into one dataframe easier

MIP_BAR_trajec <- data.frame(deSolve::ode(y = init_states,
                                          times = times,
                                          func = MIP_BAR_ode,
                                          parms = parms))

all_trajecs_df <- rbind(MIP_trajec, BAR_trajec, MIP_BAR_trajec)
all_trajecs_df$Model <- rep(c("MIP", "BAR", "MIP-BAR"), each = length(times))


#### Plot results using ggplot2 (Fig 4a in paper) ####

(time_series_matched_plot <- ggplot2::ggplot(data = all_trajecs_df, aes(x = time, y = I, col = Model)) +
   geom_line() +
   labs(x = "Time (days)", 
        y = "Proportion of infected plants (I/H)") +
   theme_bw() +
   theme(legend.position = c(0.8, 0.2)) +
   scale_color_manual(values = c("BAR" = "#F8766D", "MIP-BAR" = "black", "MIP" = "#619CFF"),
                      labels = c("BAR" = "BAR", "MIP-BAR" = "MIP-BAR, Z(0) = 0", "MIP" = "MIP"),
                      breaks = c("BAR", "MIP-BAR", "MIP"))
)
