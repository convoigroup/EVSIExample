############################################################
# SET SEED AND LOAD REQUIRED PACKAGES
############################################################

# To make sure the results can be replicated
set.seed(1234)

############################################################
# Check and load required missing packages
############################################################
if(!require(pacman)){               
  install.packages("pacman")
  library(pacman)
}

if(!require(BiocManager)){               
  install.packages("BiocManager")
  library(BiocManager)
}

p_install(BCEA, force=FALSE)             # Bayesian cost-effectiveness Analysis
p_load(BCEA)

p_install(ggplot2, force=FALSE)          # A plotting package
p_load(ggplot2)

p_install(reshape, force=FALSE)          # Data manipulation package
p_load(reshape)

p_install(plyr, force=FALSE)             # Data manipulation package
p_load(plyr)

p_install(R2jags, force=FALSE)           # Using R to run 'JAGS'
p_load(R2jags)

p_install(R2OpenBUGS, force=FALSE)       # Interfaces R with "OpenBUGS"
p_load(R2OpenBUGS)

p_install(mgcv, force=FALSE)             # Mixed GAM computation vehicle with 
p_load(mgcv)                               # automatic smoothness estimation

if(!require(EVSI)){                      # Expected value of sample information
  # To install, use devtools::install_github("annaheath/EVSI")
  devtools::install_github("annaheath/EVSI") 
  library(EVSI)
}


############################################################
# Load helper functions
############################################################

##### Generate LogNormal distribution parameters #####
lognPar <- function(m,s) {
  # m - the mean
  # s - the standard deviation
  s2 <- s^2
  mulog <- log(m) - 0.5 * log(1 + s2/m^2)
  s2log <- log(1 + (s2/m^2))
  sigmalog <- sqrt(s2log)
  # return lognormal distribution parameters
  return(params = list(mulog = mulog, sigmalog = sigmalog))
}


##### Generate Beta distribution parameters #####
betaPar <- function(mu, var) {
  # mu - the mean
  # var - the standard deviation
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  # return beta distribution parameters
  return(params = list(alpha = alpha, beta = beta))
}


############################################################
# Set literature-based parameter values
############################################################

##### Safety data (Given standard-of-care) #####
num.pat <- 111           # Number of patients (observed data)  
num.se <- 27             # Number of patients with side effects
# Number of patient with hospital care have side effect given stardard-of-care
num.home <- 10           
num.hosp <- 17           # Number of patients require hospital care
num.dead <- 1            # Number of Deaths


##### Informative prior for risk reduction #####
m.rho <- 0.65            # Mean
s.rho <- 0.1             # Standard deviation
tau.rho <- 1/s.rho^2     # Precision


##### Priors for daily probability of recovery #####

# PSA distributions for lambda.home.rec.TH (mean of beta distribution)
m.home.rec <- 0.45
# PSA distributions for lambda.home.rec.TH (variance of beta distribution)
v.home.rec <- 0.02          
# Shape parameters of beta distribution
p1.home.rec <- betaPar(m.home.rec,v.home.rec)$alpha 
# Shape parameters of beta distribution
p2.home.rec <- betaPar(m.home.rec,v.home.rec)$beta   

# PSA distributions for lambda.hosp.rec.TH (mean of beta distribution)
m.hosp.rec <- 0.35  
# PSA distributions for lambda.hosp.rec.TH (variance of beta distribution)
v.hosp.rec <- 0.02  
# Shape parameters of beta distribution
p1.hosp.rec <- betaPar(m.hosp.rec,v.hosp.rec)$alpha
# Shape parameters of beta distribution
p2.hosp.rec <- betaPar(m.hosp.rec,v.hosp.rec)$beta   


##### Prior and parameters for costs #####

# Home care 
# Mean before transformation
mu.home <- 2300   
# Standard deviation before transformation
sd.home <- 90                    
# Mean of log normal distribution
m.home <- lognPar(mu.home,sd.home)$mulog
# SD of log normal distribution
s.home <- lognPar(mu.home,sd.home)$sigmalog
# Precision of log normal distribution
tau.home <- 1/s.home^2                               

# Hospitalization
mu.hosp <- 6500
sd.hosp <- 980
m.hosp <- lognPar(mu.hosp,sd.hosp)$mulog
s.hosp <- lognPar(mu.hosp,sd.hosp)$sigmalog
tau.hosp <- 1/s.hosp^2

# Death
mu.dead <- 4200
sd.dead <- 560
m.dead <- lognPar(mu.dead,sd.dead)$mulog
s.dead <- lognPar(mu.dead,sd.dead)$sigmalog
tau.dead <- 1/s.dead^2

# Drug costs
c.drug <- c(110, 420)


##### Prior and parameters for effects/utilities #####
# QALY on chemotherapy
mu.chemo <- 0.98
var.chemo <- 0.001
p1.chemo <- betaPar(mu.chemo,var.chemo)$alpha
p2.chemo <- betaPar(mu.chemo,var.chemo)$beta

# QALY on home care
mu.e.home <- 0.5
var.e.home <- 0.02
p1.home <- betaPar(mu.e.home,var.e.home)$alpha
p2.home <- betaPar(mu.e.home,var.e.home)$beta

# QALY on hospitalization
mu.e.hosp <- 0.2
var.e.hosp <- 0.03
p1.hosp <- betaPar(mu.e.hosp,var.e.hosp)$alpha
p2.hosp <- betaPar(mu.e.hosp,var.e.hosp)$beta


##### Number of patients expected to receive treatment for predictive modeling #####
n.pred <- 1000


##### Time horizon of side effects markov model #####
TH <- 15


############################################################
# Generate PSA parameters          
############################################################

### Specify a Bayesian decision model (the Chemotherapy model) for run with JAGS
model <- function(){
  
  ########################### Side Effect Analysis ############################
  
  # Sampling distribution 
  # Number of patients who experienced side effects
  num.se ~ dbin(pi[1], num.pat)
  
  ##### Prior distribution #####
  
  # Probability of side effects given Soc
  pi[1] ~ dbeta(1, 1)                   
  # Relative risk of side effects given new treatment compared to SoC
  rho ~ dnorm(m.rho, tau.rho)
  # Reduction in probability of side effect given new treatment
  # 1-rho  
  # Probability of side effects given new treatment 
  pi[2] <- rho * pi[1]                  
  
  ##### Sampling distribution #####
  
  # Patients experienced side effects who require hospital care
  num.hosp ~ dbin(gamma.hosp, num.se)
  # Patients experienced side effects who require home care
  num.home <- num.se - num.hosp         
  
  ##### Prior distribution #####
  
  # Probability that side effects with hospital care
  gamma.hosp ~ dbeta(1, 1)
  # Probability that side effects with home care
  gamma.home <- 1 - gamma.hosp          
  
  ##### Sampling distribution #####
  
  # Patients experienced side effects who dead eventually
  num.dead ~ dbin(gamma.dead, num.hosp) 
  
  ##### Prior distribution #####
  
  # Probability that side effects result in death eventually
  gamma.dead ~ dbeta(1,4)               
  
  # Daily rate of recovery for home care and hospital care
  recover.home <- -log(1 - lambda.home.rec.TH)
  recover.hosp <- -log(1 - lambda.hosp.rec.TH)
  
  ################################ Recovery ###################################     
  
  ##### Daily probability of recovery #####
  
  # Daily probability of recovery given remain in home care
  lambda.home.rec.TH ~ dbeta(p1.home.rec, p2.home.rec)  
  # Daily probability of recovery given hospital care
  lambda.hosp.rec.TH ~ dbeta(p1.hosp.rec, p2.hosp.rec)                
  
  ##### Transition probability from home care state ##### 
  
  # Probability of staying in home care
  lambda.home.home <- 1 - lambda.home.rec - lambda.home.hosp  
  # Probability of recovering with home care
  lambda.home.rec <- (1 - lambda.home.hosp)*lambda.home.rec.TH 
  # Cumulative hazard of hospitalization
  cum.haz.hosp <- -log(1 - gamma.hosp)  
  # Daily cumulative hazard of hospitalization
  daily.cum.haz.hosp <- cum.haz.hosp/TH 
  # Daily probability of switching to hospital care
  lambda.home.hosp <- 1 - exp(-daily.cum.haz.hosp)                   
  
  ##### Transition probability from hospital care state ##### 
  
  # Probability of staying in hospital care
  lambda.hosp.hosp <- 1 - lambda.hosp.rec - lambda.hosp.dead  
  # Probability of recovering with hospital care
  lambda.hosp.rec <- (1 - lambda.hosp.dead)*lambda.hosp.rec.TH  
  # Cumulative hazard of hospitalization
  cum.haz.dead <- -log(1 - gamma.dead) 
  # Daily cumulative hazard of hospitalization
  daily.cum.haz.dead <- cum.haz.dead/TH 
  # Daily probability of death
  lambda.hosp.dead <- 1 - exp(-daily.cum.haz.dead)                   
  
  ################################# Costs #####################################
  # These are sampled direct from distributions, so no prior is needed
  c.home ~ dlnorm(m.home, tau.home)      # Cost of home care 
  c.hosp ~ dlnorm(m.hosp, tau.hosp)      # Cost of hospitalization
  c.dead ~ dlnorm(m.dead, tau.dead)      # Cost of death
  
  ################################ Effects ####################################
  # Patients who do not experience adverse events have a QoL measure of q (e.chemo) 
  e.chemo ~ dbeta(p1.chemo, p2.chemo)    # The QoL weight with no adverse events
  e.home ~ dbeta(p1.home, p2.home)       # Effect (QALY score) of home care
  e.hosp ~ dbeta(p1.hosp, p2.hosp)       # Effect (QALY score) of hospital care
  
  ############### Predictive distributions on clinical outcomes ############### 
  for (t in 1:2) {
    # Expected number of patients with side effects given Soc (SE[1]) 
    # Expected number of patients with side effects given new treatment (SE[2])
    SE[t] ~ dbin(pi[t], N)                
  }                                      
  
  ################################### VoI ####################################
  # Probability of side effects given Soc (pi1)
  pi1 <- pi[1]                           
  # Probability of side effects given new treatment (pi2)
  pi2 <- pi[2]                           
  
} 


############################################################
# Data and Parameters for the JAGS Model
############################################################  

############################ Simulation sample size ###########################
Size.Prior <- 10000

###################### Prior parameters for JAGS model ########################     
data <- list(num.pat = num.pat, 
             num.se = num.se,
             num.hosp = num.hosp, 
             num.dead = num.dead,
             p1.home.rec = p1.home.rec, 
             p2.home.rec = p2.home.rec,
             p1.hosp.rec = p1.hosp.rec, 
             p2.hosp.rec = p2.hosp.rec,
             m.rho = m.rho, 
             tau.rho = tau.rho,
             m.home = m.home, 
             tau.home = tau.home,
             m.hosp = m.hosp, 
             tau.hosp = tau.hosp,
             m.dead = m.dead, 
             tau.dead = tau.dead,
             N = n.pred,
             p1.chemo = p1.chemo, 
             p2.chemo = p2.chemo,
             p1.home = p1.home, 
             p2.home = p2.home,
             p1.hosp = p1.hosp, 
             p2.hosp = p2.hosp,
             TH = TH)

#################### The initial values generator function ####################

inits <- function(){
  list(pi = c(runif(1), NA)                       # Non-informative prior
  )
}

size.prior <- 10000                               # Size of prior
n.chains <- 3                                     # Number of chains
n.burnin <- 1000                                  # Number of burn in iterations
n.iter <- ceiling(Size.Prior/n.chains) + n.burnin # Number of iterations per chain


################# Choose the parameters in the model to monitor ###############
parameters.to.save <- c("pi1", 
                        "pi2", 
                        "rho", 
                        "gamma.hosp",
                        "gamma.dead", "SE",
                        "lambda.home.home", 
                        "lambda.home.hosp", 
                        "lambda.home.rec",
                        "lambda.hosp.hosp", 
                        "lambda.hosp.rec",
                        "lambda.hosp.dead",
                        "lambda.home.rec.TH", 
                        "lambda.hosp.rec.TH",
                        "c.home", 
                        "c.hosp", 
                        "c.dead",
                        "e.chemo", 
                        "e.home", 
                        "e.hosp",
                        "pi[1]", 
                        "pi[2]",
                        "recover.home", 
                        "recover.hosp",
                        "gamma.home")

######################## Generate Bayesian model file #########################
filein <- file.path(tempdir(), fileext="psitemp.txt")
write.model(model, filein)


##################### Run JAGS model for MCMC simulation ######################

prior.model <- jags(
  data =  data,
  inits = inits,
  parameters.to.save = parameters.to.save,
  model.file = filein, 
  n.chains = n.chains, 
  n.iter = n.iter, 
  n.thin = 1, 
  n.burnin = n.burnin,
  progress.bar = "none"
) 



############################################################
# Run the Markov model simulation
############################################################

# Calculate average time each patient with adverse events spends in the health 
  #states of the Markov model
markov.model <- function(SE,
                         lambda.home.home, 
                         lambda.home.hosp,
                         lambda.home.rec,
                         lambda.hosp.hosp, 
                         lambda.hosp.rec, 
                         lambda.hosp.dead, 
                         TH)
  # SE: predictive distribution of number of patients in 1000 that would 
        #experienced adverse events, which is a two-vector containing the number
        #of adverse events for Soc and novel treatment
  
{ # MM.mat: Construct a transition matrix that governs the states transition 
  # Row of MM.mat: Four Markov states (Home care, Hospital care, Recovery, Death)
  # Column of MM.mat: transition probability from one state to another
  MM.mat <- matrix(c(lambda.home.home, lambda.home.hosp, lambda.home.rec, 0,
                     0, lambda.hosp.hosp, lambda.hosp.rec, lambda.hosp.dead,
                     0, 0, 1, 0,
                     0, 0, 0, 1),
                   nrow = 4,
                   ncol = 4,
                   byrow = TRUE)
  
  # Trace matrix: calculate number of patients in each state for each time point 
  # It has 3 dimensions: # number of states,
                         # number of time points,
                         # length of side effects
  trace <- array(0, dim = c(4, TH + 1, length(SE))) 
  trace[1, 1, ] <- SE      
  
  for(i in 2:(TH + 1)){   
    trace[, i, 1] <- trace[, i - 1, 1] %*% MM.mat 
    trace[, i, 2] <- trace[, i - 1, 2] %*% MM.mat
  }
  
  return(trace)    # return "trace" which is a 4*16*2 array
}


############################################################
# Estimate cost
############################################################
########################### Costs for two treatments ##########################
# Calculated based on costs for each health state
costs <- function(SE,
                  lambda.home.home, 
                  lambda.home.hosp, 
                  lambda.home.rec,
                  lambda.hosp.hosp, 
                  lambda.hosp.rec, 
                  lambda.hosp.dead,
                  c.home, 
                  c.hosp, 
                  c.dead,
                  N = n.pred, 
                  time.horiz = TH){
  trace <- markov.model(SE,
                        lambda.home.home, 
                        lambda.home.hosp, 
                        lambda.home.rec,
                        lambda.hosp.hosp, 
                        lambda.hosp.rec, 
                        lambda.hosp.dead, 
                        time.horiz)
  
  # Costs of four states
  c.states <- c(c.home, c.hosp, 0, 0)
  
  # Cost of Side Effects (for Soc and novel treatment) 
  c.se <- array(NA, dim = 2)
  
  # Average cost for both Soc and novel treatment per person per time point
  # This includes one-off cost of death for patients who died at the end of 15 days
  c.se[1] <- (sum(c.states %*% trace[,,1]) + c.dead * trace[4, TH + 1, 1])/
    (N * (TH + 1)) 
  c.se[2] <- (sum(c.states %*% trace[,,2]) + c.dead * trace[4, TH + 1, 2])/
    (N * (TH + 1)) 
  
  # Cost for both Soc and novel treatment   
  c.drug <- c(110, 420)
  cost <- c(c.drug + c.se)
  return(cost)
}


##### Baseline PSA simulations results for costs #####
c <- matrix(NA, ncol = 2, nrow = size.prior)
for(l in 1:size.prior){
  c[l,] <- costs(prior.model$BUGSoutput$sims.list[["SE"]][l, ],
                 prior.model$BUGSoutput$sims.list[["lambda.home.home"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.home.hosp"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.home.rec"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.hosp.hosp"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.hosp.rec"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.hosp.dead"]][l],
                 prior.model$BUGSoutput$sims.list[["c.home"]][l],
                 prior.model$BUGSoutput$sims.list[["c.hosp"]][l],
                 prior.model$BUGSoutput$sims.list[["c.dead"]][l])
}



############################################################
# Estimate utility
############################################################

######################## Effectiveness for two treatments #####################
# Calculated based on QoL measure for each health state
effects <- function(SE,
                    lambda.home.home, 
                    lambda.home.hosp, 
                    lambda.home.rec,
                    lambda.hosp.hosp, 
                    lambda.hosp.rec, 
                    lambda.hosp.dead,
                    e.chemo, 
                    e.home, 
                    e.hosp,
                    N = n.pred ,
                    time.horiz = TH){
  trace <- markov.model(SE,
                        lambda.home.home, 
                        lambda.home.hosp, 
                        lambda.home.rec,
                        lambda.hosp.hosp, 
                        lambda.hosp.rec, 
                        lambda.hosp.dead, 
                        time.horiz)
  
  # Effectiveness for four states
  e.states <- c(e.home, e.hosp, e.chemo, 0)
  
  # Total QALY of side effects for both Soc and novel treatment
  q.se <- array(NA, dim = 2)
  q.se[1] <- sum(e.states %*% trace[,,1])
  q.se[2] <- sum(e.states %*% trace[,,2])
  
  # QALY of total number of patients who do not experience adverse events for 15 days
  q.chemo <- (N - SE) * e.chemo * (TH + 1)
  
  # Average effect for both Soc and novel treatment per person per time point
  effs<-c(q.se + q.chemo)/
    (N * (TH + 1))           
  return(effs)
}



##### Baseline PSA simulations results for utility (effectivenes) #####
e <- matrix(NA, ncol = 2, nrow = size.prior)
for(l in 1:size.prior){
  e[l,] <- effects(prior.model$BUGSoutput$sims.list[["SE"]][l, ], 
                   prior.model$BUGSoutput$sims.list[["lambda.home.home"]][l],
                   prior.model$BUGSoutput$sims.list[["lambda.home.hosp"]][l],
                   prior.model$BUGSoutput$sims.list[["lambda.home.rec"]][l],
                   prior.model$BUGSoutput$sims.list[["lambda.hosp.hosp"]][l],
                   prior.model$BUGSoutput$sims.list[["lambda.hosp.rec"]][l],
                   prior.model$BUGSoutput$sims.list[["lambda.hosp.dead"]][l],
                   prior.model$BUGSoutput$sims.list[["e.chemo"]][l],
                   prior.model$BUGSoutput$sims.list[["e.home"]][l],
                   prior.model$BUGSoutput$sims.list[["e.hosp"]][l])
}


############################################################
# Perform baseline cost and effectiveness analysis
############################################################
mod <- bcea(e, c, ref = 2, interventions = c("Standard of Care", 
                                             "New Treatment" ))

# Produces a scatter plot of the cost-effectiveness plane
plot.stand <- ceplane.plot(mod, wtp = 30000, graph = "ggplot")
plot.stand[[2]] <- plot.stand[[2]][-6]
plot.stand

# Performs baseline cost-effectiveness analysis for specified willingness to pay
wtp <- 30000
NB <- e*wtp - c


############################################################
# Package cost, utility and model parameters for output
############################################################
# Matrix of parameters of interest from baseline model
extra.lines <- (Size.Prior + 1):dim(prior.model$BUGSoutput$sims.matrix)[1]
theta <- as.data.frame(prior.model$BUGSoutput$sims.matrix[-extra.lines,c("pi1",
                                                                         "rho",
                                                                         "gamma.hosp",
                                                                         "gamma.dead",
                                                                         "lambda.home.rec.TH",
                                                                         "lambda.hosp.rec.TH",
                                                                         "SE[1]",
                                                                         "SE[2]",
                                                                         "recover.home",
                                                                         "recover.hosp",
                                                                         "pi2",
                                                                         "gamma.home")])
colnames(theta) <- c("pi1",
                     "rho",
                     "gamma.hosp",
                     "gamma.dead",
                     "lambda.home.rec.TH",
                     "lambda.hosp.rec.TH",
                     "SE[1]",
                     "SE[2]",
                     "recover.home",
                     "recover.hosp",
                     "pi2",
                     "gamma.home")


############################################################
# Expected value of perfect information
############################################################
# Calculate EVPI for a specified willingness to pay
mod$evi[which(mod$k == 30000)] 

# From Nested Simulations avaliable on request
# evsi.true <- 21.10214  

