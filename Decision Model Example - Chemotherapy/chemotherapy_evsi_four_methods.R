###################################################################################
# Chemotherapy
###################################################################################

############################################################
# Load model and EVSI methods scripts
############################################################
source(".//Decision Model Example - Chemotherapy//chemotherapy_decision_model.R")
source(".//EVSI Estimation Methods//evsi_moment_matching_method.R")
source(".//EVSI Estimation Methods//evsi_Gaussian_approximation_method.R")
source(".//EVSI Estimation Methods//evsi_importance_sampling_method.R")
source(".//EVSI Estimation Methods//evsi_regression_based_method.R")

############################################################
# Moment Matching Method
############################################################

decision_model_cea <- list(cost=c,
                           effect=e,
                           cost_funtion=costs,
                           effect_funtion=effects,
                           wtp=30000)

decision_model_psa <- theta[,c("pi1", 
                               "rho",
                               "gamma.hosp", 
                               "gamma.dead",
                               "lambda.home.rec.TH", 
                               "lambda.hosp.rec.TH",
                               "recover.home", 
                               "recover.hosp",
                               "pi2", 
                               "gamma.home")]

model.data.full <- function(){
  # Number of side effects in standard care
  X.SE1 ~ dbin(pi[1], n) 
  # Number of side effects in novel treatment
  X.SE2 ~ dbin(pi[2], n)  
  # Total number of patients in hospital care
  X.N.hosp ~ dbinom(gamma.hosp, X.SE1 + X.SE2)
  # Total number of patients dead
  X.N.dead ~ dbin(gamma.dead, X.N.hosp)        
  
  # Daily rate of recovery for home care and hospital care
  recover.home <- -log(1-lambda.home.rec.TH)
  recover.hosp <- -log(1-lambda.hosp.rec.TH)
  
  # Recovery time at home
  for(i in 1:N.home){
    T.rec.home[i] ~ dexp(recover.home)
  }
  
  # Recovery time at hospital
  for(i in 1:N.hosp){
    T.rec.hosp[i] ~ dexp(recover.hosp)
  }
  
  # Number of side effects
  num.se ~ dbin(pi[1], num.pat)    
  pi[1] ~ dbeta(1, 1)               
  
  rho ~ dnorm(m.rho, tau.rho)       
  pi[2] <- rho*pi[1]     
  
  for (t in 1:2) {
    SE[t] ~ dbin(pi[t], N)
  }
  
  ## Number of hopistalized individual (Sampling distribution)
  
  # Patients experienced side effects who require hospital care
  num.hosp ~ dbin(gamma.hosp, num.se)   
  # Probability that side effects with hospital care
  gamma.hosp ~ dbeta(1, 1)                    
  
  ## Number of death
  
  # Patients experienced side effects who dead eventually
  num.dead ~ dbin(gamma.dead, num.hosp)
  # Probability that side effects result in death eventually 
  gamma.dead ~ dbeta(1, 4)                     
  
  ## Daily probability of recovery
  lambda.home.rec.TH ~ dbeta(p1.home.rec, p2.home.rec)
  lambda.hosp.rec.TH ~ dbeta(p1.hosp.rec, p2.hosp.rec)
  
  ## Transition probability from home care state
  
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
  
  ## Transition probability from hospital care state
  
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
  
  # Costs
  c.home ~ dlnorm(m.home, tau.home)  # Cost of home care 
  c.hosp ~ dlnorm(m.hosp, tau.hosp)  # Cost of hospitalization
  c.dead ~ dlnorm(m.dead, tau.dead)  # Cost of death
  
  # Effects
  e.chemo ~ dbeta(p1.chemo,p2.chemo) # The QoL weight with no adverse events
  e.home ~ dbeta(p1.home,p2.home)    # Effect (QALY score) of home care
  e.hosp ~ dbeta(p1.hosp,p2.hosp)    # Effect (QALY score) of hospital care 
  
  # Probability of side effects given Soc (pi1)
  pi1 <- pi[1] 
  # Probability of side effects given new treatment (pi2)
  pi2 <- pi[2]                          
}
inits <- function(){
  list(pi = c(runif(1), NA))                       # Non-informative prior
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=data,
                            update_prior_parameter=c("X.SE1", 
                                                     "X.SE2",
                                                     "X.home", 
                                                     "X.hosp",
                                                     "X.N.dead", 
                                                     "X.N.hosp",
                                                     "N.home", 
                                                     "N.hosp",
                                                     "T.rec.home", 
                                                     "T.rec.hosp"),
                            inits=inits,
                            iter_number=10000,
                            burnin_number=3000,
                            mcmc_parameter=c("pi1",
                                             "rho",
                                             "gamma.hosp",
                                             "gamma.dead", 
                                             "lambda.home.rec.TH", 
                                             "lambda.hosp.rec.TH"))

data_generation_function <- function(n, theta){
  X.SE1 <- rbinom(1, n, theta$pi1)
  X.SE2 <- rbinom(1, n, theta$pi2)
  
  X.N.hosp <- rbinom(1,X.SE1 + X.SE2, theta$gamma.hosp)
  X.N.dead <- rbinom(1,X.N.hosp, theta$gamma.dead)
  
  N.home <- X.SE1 + X.SE2 - X.N.hosp
  T.rec.home <- rexp(N.home, theta$recover.home)
  X.home <- sum(T.rec.home)
  
  N.hosp <- X.N.hosp - X.N.dead
  T.rec.hosp <- rexp(N.hosp, theta$recover.hosp)
  X.hosp <- sum(T.rec.hosp)
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list(X.SE1=X.SE1, 
                           X.SE2=X.SE2,
                           X.home=X.home, 
                           X.hosp=X.hosp,
                           X.N.dead=X.N.dead, 
                           X.N.hosp=X.N.hosp,
                           N.home=N.home, 
                           N.hosp=N.hosp,
                           T.rec.home=T.rec.home, 
                           T.rec.hosp=T.rec.hosp)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  
  return(D)
}
evsi_specification <- list(evppi_parameter=c("pi1",
                                             "rho",
                                             "gamma.hosp",
                                             "gamma.dead", 
                                             "lambda.home.rec.TH",
                                             "lambda.hosp.rec.TH"),
                           nested_sample_size=50,
                           quantile_parameter=c("pi1", 
                                                "rho",
                                                "gamma.hosp", 
                                                "gamma.dead",
                                                "lambda.home.rec.TH", 
                                                "lambda.hosp.rec.TH",
                                                "recover.home", 
                                                "recover.hosp",
                                                "pi2", 
                                                "gamma.home"),
                           sim_data_size=c(150))

evsi_moment_matching(decision_model_cea,
                     decision_model_psa,
                     decision_model_mcmc,
                     data_generation_function,
                     evsi_specification)

############################################################
# Gaussian Approximation Method
############################################################

decision_model_cea <- list(cost=c,
                           effect=e,
                           wtp=30000)

decision_model_psa <- theta[,c("pi1", 
                               "rho",
                               "gamma.hosp", 
                               "gamma.dead",
                               "lambda.home.rec.TH", 
                               "lambda.hosp.rec.TH",
                               "recover.home", 
                               "recover.hosp",
                               "pi2", 
                               "gamma.home")]
model.data.full <- function(){
  # Number of side effects in standard care
  X.SE1 ~ dbin(pi[1], n)       
  # Number of side effects in novel treatment
  X.SE2 ~ dbin(pi[2], n)               
  # Total number of patients in hospital care
  X.N.hosp ~ dbinom(gamma.hosp, X.SE1 + X.SE2) 
  # Total number of patients dead
  X.N.dead ~ dbin(gamma.dead, X.N.hosp)        
  
  # Daily rate of recovery for home care and hospital care
  recover.home <- -log(1-lambda.home.rec.TH)
  recover.hosp <- -log(1-lambda.hosp.rec.TH)
  
  # Recovery time at home
  for(i in 1:N.home){
    T.rec.home[i] ~ dexp(recover.home)
  }
  
  # Recovery time at hospital
  for(i in 1:N.hosp){
    T.rec.hosp[i] ~ dexp(recover.hosp)
  }
  
  # Number of side effects
  num.se ~ dbin(pi[1], num.pat)    
  pi[1] ~ dbeta(1, 1)               
  
  rho ~ dnorm(m.rho, tau.rho)       
  pi[2] <- rho*pi[1]     
  
  for (t in 1:2) {
    SE[t] ~ dbin(pi[t], N)
  }
  
  ## Number of hopistalized individual (Sampling distribution)
  # Patients experienced side effects who require hospital care
  num.hosp ~ dbin(gamma.hosp, num.se)    
  # Probability that side effects with hospital care
  gamma.hosp ~ dbeta(1, 1)                    
  
  ## Number of death
  # Patients experienced side effects who dead eventually
  num.dead ~ dbin(gamma.dead, num.hosp)  
  # Probability that side effects result in death eventually
  gamma.dead ~ dbeta(1, 4)                      
  
  ## Daily probability of recovery
  lambda.home.rec.TH ~ dbeta(p1.home.rec, p2.home.rec)
  lambda.hosp.rec.TH ~ dbeta(p1.hosp.rec, p2.hosp.rec)
  
  ## Transition probability from home care state
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
  
  ## Transition probability from hospital care state
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
  
  ## Costs
  # Cost of home care 
  c.home ~ dlnorm(m.home, tau.home)
  # Cost of hospitalization
  c.hosp ~ dlnorm(m.hosp, tau.hosp) 
  # Cost of death
  c.dead ~ dlnorm(m.dead, tau.dead)      
  
  ## Effects
  # The QoL weight with no adverse events
  e.chemo ~ dbeta(p1.chemo,p2.chemo)
  # Effect (QALY score) of home care
  e.home ~ dbeta(p1.home,p2.home)   
  # Effect (QALY score) of hospital care
  e.hosp ~ dbeta(p1.hosp,p2.hosp)         
  
  # Probability of side effects given Soc (pi1)
  pi1 <- pi[1]                  
  # Probability of side effects given new treatment (pi2)
  pi2 <- pi[2]                          
}
inits <- function(){
  list(pi = c(runif(1), NA))                  # Non-informative prior
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=data,
                            update_prior_parameter=c("X.SE1", 
                                                     "X.SE2",
                                                     "X.home", 
                                                     "X.hosp",
                                                     "X.N.dead", 
                                                     "X.N.hosp",
                                                     "N.home", 
                                                     "N.hosp",
                                                     "T.rec.home", 
                                                     "T.rec.hosp"),                            
                            inits=inits,
                            chain_number=3,
                            burnin_number=1000,
                            mcmc_parameter=c("pi1","rho",
                                             "gamma.hosp","gamma.dead",
                                             "lambda.home.rec.TH","lambda.hosp.rec.TH"))

data_generation_function <- function(n, theta){
  X.SE1 <- rbinom(1, n, theta$pi1)
  X.SE2 <- rbinom(1, n, theta$pi2)
  
  X.N.hosp <- rbinom(1,X.SE1 + X.SE2, theta$gamma.hosp)
  X.N.dead <- rbinom(1,X.N.hosp, theta$gamma.dead)
  
  N.home <- X.SE1 + X.SE2 - X.N.hosp
  T.rec.home <- rexp(N.home, theta$recover.home)
  X.home <- sum(T.rec.home)
  
  N.hosp <- X.N.hosp - X.N.dead
  T.rec.hosp <- rexp(N.hosp, theta$recover.hosp)
  X.hosp <- sum(T.rec.hosp)
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list(X.SE1=X.SE1, X.SE2=X.SE2,
                           X.home=X.home, X.hosp=X.hosp,
                           X.N.dead=X.N.dead, X.N.hosp=X.N.hosp,
                           N.home=N.home, N.hosp=N.hosp,
                           T.rec.home=T.rec.home, T.rec.hosp=T.rec.hosp)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}
evsi_specification <- list(evppi_parameter=c("pi1","rho",
                                             "gamma.hosp","gamma.dead",
                                             "lambda.home.rec.TH","lambda.hosp.rec.TH"),
                           economic_model_size=30,
                           outer_cond_exp_loop_size=1000,
                           inner_cond_exp_loop_size=10000,
                           sim_data_size=c(150),
                           mcmc_retry=5)

evsi_gaussian_approximation(decision_model_cea,
                            decision_model_psa,
                            decision_model_mcmc,
                            data_generation_function,
                            evsi_specification)

############################################################
# Importance Sampling Method
############################################################

decision_model_cea <- list(cost=c,
                           effect=e,
                           wtp=30000)

decision_model_psa <- theta[,c("pi1", 
                               "rho",
                               "gamma.hosp", 
                               "gamma.dead",
                               "lambda.home.rec.TH", 
                               "lambda.hosp.rec.TH",
                               "recover.home", 
                               "recover.hosp",
                               "pi2", 
                               "gamma.home")]

data_generation_function <- function(n, theta){
  X.SE1 <- rbinom(1, n, theta$pi1)
  X.SE2 <- rbinom(1, n, theta$pi2)
  
  X.N.hosp <- rbinom(1,X.SE1 + X.SE2, theta$gamma.hosp)
  X.N.dead <- rbinom(1,X.N.hosp, theta$gamma.dead)
  
  N.home <- X.SE1 + X.SE2 - X.N.hosp
  T.rec.home <- rexp(N.home, theta$recover.home)
  X.home <- sum(T.rec.home)
  
  N.hosp <- X.N.hosp - X.N.dead
  T.rec.hosp <- rexp(N.hosp, theta$recover.hosp)
  X.hosp <- sum(T.rec.hosp)
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list(X.SE1=X.SE1, 
                           X.SE2=X.SE2,
                           X.home=X.home, 
                           X.hosp=X.hosp,
                           X.N.dead=X.N.dead, 
                           X.N.hosp=X.N.hosp,
                           N.home=N.home, 
                           N.hosp=N.hosp,
                           T.rec.home=T.rec.home, 
                           T.rec.hosp=T.rec.hosp)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

likelihood_compute_function <- function(n, D, theta){
  PSA.sim <- dim(theta)[1]
  colnames(D) <- c("SE1", 
                   "SE2", 
                   "home", 
                   "hosp", 
                   "N.dead", 
                   "N.hosp")
  colnames(theta) <- c("pi1", 
                       "rho", 
                       "gamma.hosp", 
                       "gamma.dead",
                       "lambda.home.rec.TH", 
                       "lambda.hosp.rec.TH",
                       "home", 
                       "hosp",
                       "pi2", 
                       "gamma.home")
  
  l <- dbinom(D[,"SE1"], n, theta[,"pi1"], log = TRUE) +
    dbinom(D[,"SE2"], n, theta[,"pi2"], log = TRUE) +
    dbinom(D[,"N.hosp"], D[,"SE1"] + D[,"SE2"], theta[,"gamma.hosp"], log = TRUE) +
    dbinom(D[,"N.dead"], D[,"N.hosp"], theta[,"gamma.dead"], log = TRUE) +
    log(theta[,"home"]) * (D[,"SE1"] + D[,"SE2"] - D[,"N.hosp"]) - theta[,"home"] * D[,"home"] +
    log(theta[,"hosp"]) * (D[,"N.hosp"] - D[,"N.dead"]) - theta[,"hosp"] * D[,"hosp"]
  
  return(exp(l))
}

evsi_specification <- list(evppi_parameters=c("pi1",
                                              "rho",
                                              "gamma.hosp",
                                              "gamma.dead",
                                              "lambda.home.rec.TH",
                                              "lambda.hosp.rec.TH"),
                           gam_tensor_prod_dim=3,
                           likelihood_data_variables=c("X.SE1", 
                                                       "X.SE2", 
                                                       "X.home", 
                                                       "X.hosp", 
                                                       "X.N.dead", 
                                                       "X.N.hosp"),
                           outer_cond_exp_loop_size=10000,
                           sim_data_size=c(150))

evsi_importance_sample(decision_model_cea,
                       decision_model_psa,
                       data_generation_function,
                       likelihood_compute_function,
                       evsi_specification)

############################################################
# Regression-Based Method
############################################################

decision_model_cea <- list(cost=c,
                           effect=e,
                           wtp=30000)

decision_model_psa <- theta[,c("pi1", 
                               "rho",
                               "gamma.hosp", 
                               "gamma.dead",
                               "lambda.home.rec.TH", 
                               "lambda.hosp.rec.TH",
                               "recover.home", 
                               "recover.hosp",
                               "pi2", 
                               "gamma.home")]

data_generation_function <- function(n, theta){
  X.SE1 <- rbinom(1, n, theta$pi1)
  X.SE2 <- rbinom(1, n, theta$pi2)
  
  X.N.hosp <- rbinom(1,X.SE1 + X.SE2, theta$gamma.hosp)
  X.N.dead <- rbinom(1,X.N.hosp, theta$gamma.dead)
  
  N.home <- X.SE1 + X.SE2 - X.N.hosp
  T.rec.home <- rexp(N.home, theta$recover.home)
  X.home <- sum(T.rec.home)
  
  N.hosp <- X.N.hosp - X.N.dead
  T.rec.hosp <- rexp(N.hosp, theta$recover.hosp)
  X.hosp <- sum(T.rec.hosp)
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list(X.SE1=X.SE1, 
                           X.SE2=X.SE2,
                           X.home=X.home, 
                           X.hosp=X.hosp,
                           X.N.dead=X.N.dead, 
                           X.N.hosp=X.N.hosp,
                           N.home=N.home, 
                           N.hosp=N.hosp,
                           T.rec.home=T.rec.home, 
                           T.rec.hosp=T.rec.hosp)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(gam_parameters=c("X.home", 
                                            "X.SE1", 
                                            "X.SE2", 
                                            "X.N.hosp", 
                                            "X.hosp", 
                                            "X.N.dead"),
                           gam_tensor_prod_dim=3,
                           outer_cond_exp_loop_size=10000,
                           sim_data_size=c(150))

evsi_regression_based(decision_model_cea,
                      decision_model_psa,
                      data_generation_function,
                      evsi_specification)
  