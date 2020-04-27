###################################################################################
# Chronic Pain
###################################################################################

############################################################
# Load model and EVSI methods scripts
############################################################
source(".//Decision Model Example - Chronic Pain//chronic_pain_decision_model.R")
source(".//EVSI Estimation Methods//evsi_moment_matching_method.R")
source(".//EVSI Estimation Methods//evsi_Gaussian_approximation_method.R")
source(".//EVSI Estimation Methods//evsi_importance_sampling_method.R")
source(".//EVSI Estimation Methods//evsi_regression_based_method.R")

############################################################
# Moment Matching Method
############################################################

cost_funtion <- function(c.t1=costs.t1,c.t2=costs.t2){
  return(c(c.t1[i],c.t2[i]))
}

effect_funtion <- function(u.l1.noae,u.l1.withdraw.noae){
  u.l3<-(u.l1.noae+u.l1.ae[i])/2
  
  #Discontinuing treatment
  u.dist<-u.l1.withdraw.noae*0.8
  
  u.Mat<-matrix(c(u.l1.noae,
                  u.l1.ae[i],
                  u.l1.withdraw.ae[i],
                  u.l1.withdraw.noae,
                  u.l1.noae*u.l2,
                  u.l1.ae[i]*u.l2,
                  u.l1.withdraw.ae[i]*u.l2,
                  u.l1.withdraw.noae*u.l2,
                  u.l3,
                  u.dist)*7/365.25,nrow=10)
  
  effects.t1<-sum(Prob.Array.1[,,i]%*%u.Mat)
  
  effects.t2<-sum(Prob.Array.2[,,i]%*%u.Mat)
  
  return(c(effects.t1,effects.t2))
}

decision_model_cea <- list(cost=discount.15*cbind(costs.t1,costs.t2),
                           effect=discount.15*cbind(utility.t1,utility.t2),
                           cost_funtion=cost_funtion,
                           effect_funtion=effect_funtion,
                           wtp=20000)

decision_model_psa <- as.data.frame(pars)

model.data.full <- function(){
  for(i in 1:n){
    utility.noae[i]~dbeta(s.X2,r.X2)
    utility.ae[i]~dbeta(s.X4,r.X4)
  }
  
  s.X2 <- u.l1.noae*( (u.l1.noae*(1-u.l1.noae)/(sig.X.noae*sig.X.noae))-1)
  r.X2 <- (1-u.l1.noae)*((u.l1.noae*(1-u.l1.noae)/(sig.X.noae*sig.X.noae))-1)
  s.X4 <- u.l1.withdraw.noae*( (u.l1.withdraw.noae*(1-u.l1.withdraw.noae)/(sig.X.with*sig.X.with))-1)
  r.X4 <- (1-u.l1.withdraw.noae)*((u.l1.withdraw.noae*(1-u.l1.withdraw.noae)/(sig.X.with*sig.X.with))-1)
  
  #Utilities
  u.l1.noae~dbeta(s.u.l1.noae, r.u.l1.noae)
  #Withdraw due to other reasons
  u.l1.withdraw.noae~dbeta(s.u.l1.withdraw.noae, r.u.l1.withdraw.noae)
}

inits <- function(){
  list(u.l1.noae=rep(0.8,1),
       u.l1.withdraw.noae=rep(0.6,1))
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=list(),
                            update_prior_parameter=c("r.u.l1.noae", 
                                                     "r.u.l1.withdraw.noae",
                                                     "s.u.l1.noae", 
                                                     "s.u.l1.withdraw.noae",
                                                     "utility.noae", 
                                                     "utility.ae",
                                                     "sig.X.noae", 
                                                     "sig.X.with"),
                            inits=inits,
                            iter_number=10000,
                            burnin_number=1000,
                            mcmc_parameter=c("u.l1.noae", "u.l1.withdraw.noae"))

data_generation_function <- function(n, theta){
  #Utilities
  s.u.l1.noae <- beta.par(0.695000000)$alpha
  r.u.l1.noae <- beta.par(0.695000000)$beta
  
  #Withdraw due to other reasons
  s.u.l1.withdraw.noae <- beta.par(0.405000000)$alpha
  r.u.l1.withdraw.noae <- beta.par(0.405000000)$beta
  
  # The standard deviation for the data without adverse event
  sig.X.noae <- 0.300
  
  # The standard deviation for withdrawn
  sig.X.with <- 0.310
  
  # Random generation for the Beta distribution with parameters for utility of no adverse effect
  X2 <- rbeta(n,
              betaPar(theta$u.l1.noae, sig.X.noae)$a,
              betaPar(theta$u.l1.noae, sig.X.noae)$b)
  
  # Random generation for the Beta distribution with parameters for utility of adverse effect
  X4 <- rbeta(n,
              betaPar(theta$u.l1.withdraw.noae, sig.X.with)$a,
              betaPar(theta$u.l1.withdraw.noae, sig.X.with)$b)
  
  # Random generation for the binomial distribution with parameters size and probability
  missingness <- rbinom(n, 1, 0.687)
  
  X2[which(X2*missingness==0)] <- NA
  X4[which(X4*missingness==0)] <- NA
  
  # Summarize simulated data as geometric mean
  X.gam <- c(psych::geometric.mean(X2),
             psych::geometric.mean(1-X2))
  X.gam.with <- c(psych::geometric.mean(X4),
                  psych::geometric.mean(1-X4))
  
  # Output simulated data and its summary statistics
  r <- list()
  
  # Sum(rand) x 2 matrices for Beta distribution parameters 
  r$simulated_data <- list(r.u.l1.noae, 
                           r.u.l1.withdraw.noae,
                           s.u.l1.noae, 
                           s.u.l1.withdraw.noae,
                           X2, 
                           X4,
                           sig.X.noae, 
                           sig.X.with)
  
  names(r$simulated_data) <- c("r.u.l1.noae", 
                               "r.u.l1.withdraw.noae",
                               "s.u.l1.noae", 
                               "s.u.l1.withdraw.noae",
                               "utility.noae", 
                               "utility.ae",
                               "sig.X.noae", 
                               "sig.X.with")
  
  # Sum(rand) x 2 matrices for summarized statistics
  r$summary_statistics <- list("X.gam.1"=X.gam[1], 
                               "X.gam.2"=X.gam[2],
                               "X.gam.with.1"=X.gam.with[1], 
                               "X.gam.with.2"=X.gam.with[2])
  
  return(r)
}
evsi_specification <- list(evppi_parameters=c("u.l1.noae", 
                                              "u.l1.withdraw.noae"),
                           nested_sample_size=50,
                           quantile_parameter=c("u.l1.noae", 
                                                "u.l1.withdraw.noae"),
                           sim_data_size=c(150))

evsi_moment_matching(decision_model_cea,
                     decision_model_psa,
                     decision_model_mcmc,
                     data_generation_function,
                     evsi_specification)

############################################################
# Gaussian Approximation Method
############################################################

decision_model_cea <- list(cost=discount.15*cbind(costs.t1,costs.t2),
                           effect=discount.15*cbind(utility.t1,utility.t2),
                           wtp=20000)

decision_model_psa <- as.data.frame(pars)

model.data.full <- function(){
  for(i in 1:n){
    utility.noae[i]~dbeta(s.X2,r.X2)
    utility.ae[i]~dbeta(s.X4,r.X4)
  }
  
  s.X2 <- u.l1.noae*( (u.l1.noae*(1-u.l1.noae)/(sig.X.noae*sig.X.noae))-1)
  r.X2 <- (1-u.l1.noae)*((u.l1.noae*(1-u.l1.noae)/(sig.X.noae*sig.X.noae))-1)
  s.X4 <- u.l1.withdraw.noae*( (u.l1.withdraw.noae*(1-u.l1.withdraw.noae)/(sig.X.with*sig.X.with))-1)
  r.X4 <- (1-u.l1.withdraw.noae)*((u.l1.withdraw.noae*(1-u.l1.withdraw.noae)/(sig.X.with*sig.X.with))-1)
  
  #Utilities
  u.l1.noae~dbeta(s.u.l1.noae, r.u.l1.noae)
  #Withdraw due to other reasons
  u.l1.withdraw.noae~dbeta(s.u.l1.withdraw.noae, r.u.l1.withdraw.noae)
}

inits <- function(){
  list(u.l1.noae=rep(0.8,1),
       u.l1.withdraw.noae=rep(0.6,1))
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=list(),
                            update_prior_parameter=c("r.u.l1.noae", 
                                                     "r.u.l1.withdraw.noae",
                                                     "s.u.l1.noae", 
                                                     "s.u.l1.withdraw.noae",
                                                     "utility.noae", 
                                                     "utility.ae",
                                                     "sig.X.noae", 
                                                     "sig.X.with"),
                            inits=inits,
                            chain_number=1,
                            burnin_number=100,
                            mcmc_parameter=c("u.l1.noae", 
                                             "u.l1.withdraw.noae"))

data_generation_function <- function(n, theta){
  #Utilities
  s.u.l1.noae <- beta.par(0.695000000)$alpha
  r.u.l1.noae <- beta.par(0.695000000)$beta
  
  #Withdraw due to other reasons
  s.u.l1.withdraw.noae <- beta.par(0.405000000)$alpha
  r.u.l1.withdraw.noae <- beta.par(0.405000000)$beta
  
  # The standard deviation for the data without adverse event
  sig.X.noae <- 0.300
  
  # The standard deviation for withdrawn
  sig.X.with <- 0.310
  
  # Random generation for the Beta distribution with parameters for utility of no adverse effect
  X2 <- rbeta(n,
              betaPar(theta$u.l1.noae, sig.X.noae)$a,
              betaPar(theta$u.l1.noae, sig.X.noae)$b)
  
  # Random generation for the Beta distribution with parameters for utility of adverse effect
  X4 <- rbeta(n,
              betaPar(theta$u.l1.withdraw.noae, sig.X.with)$a,
              betaPar(theta$u.l1.withdraw.noae, sig.X.with)$b)
  
  # Random generation for the binomial distribution with parameters size and probability
  missingness <- rbinom(n, 1, 0.687)
  
  X2[which(X2*missingness==0)] <- NA
  X4[which(X4*missingness==0)] <- NA
  
  # Summarize simulated data as geometric mean
  X.gam <- c(psych::geometric.mean(X2),
             psych::geometric.mean(1-X2))
  X.gam.with <- c(psych::geometric.mean(X4),
                  psych::geometric.mean(1-X4))
  
  # Output simulated data and its summary statistics
  r <- list()
  
  # Sum(rand) x 2 matrices for Beta distribution parameters 
  r$simulated_data <- list(r.u.l1.noae, 
                           r.u.l1.withdraw.noae,
                           s.u.l1.noae, 
                           s.u.l1.withdraw.noae,
                           X2, 
                           X4,
                           sig.X.noae, 
                           sig.X.with)
  
  names(r$simulated_data) <- c("r.u.l1.noae", 
                               "r.u.l1.withdraw.noae",
                               "s.u.l1.noae", 
                               "s.u.l1.withdraw.noae",
                               "utility.noae", 
                               "utility.ae",
                               "sig.X.noae", 
                               "sig.X.with")
  
  # Sum(rand) x 2 matrices for summarized statistics
  r$summary_statistics <- list("X.gam.1"=X.gam[1], 
                               "X.gam.2"=X.gam[2],
                               "X.gam.with.1"=X.gam.with[1], 
                               "X.gam.with.2"=X.gam.with[2])
  
  return(r)
}

evsi_specification <- list(evppi_parameters=c("u.l1.noae", 
                                              "u.l1.withdraw.noae"),
                           economic_model_size=40,
                           outer_cond_exp_loop_size=10000,
                           inner_cond_exp_loop_size=5000,
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

decision_model_cea <- list(cost=discount.15*cbind(costs.t1,costs.t2),
                           effect=discount.15*cbind(utility.t1,utility.t2),
                           wtp=20000)

decision_model_psa <- as.data.frame(pars)

data_generation_function <- function(n, theta){
  #Utilities
  s.u.l1.noae <- beta.par(0.695000000)$alpha
  r.u.l1.noae <- beta.par(0.695000000)$beta
  
  #Withdraw due to other reasons
  s.u.l1.withdraw.noae <- beta.par(0.405000000)$alpha
  r.u.l1.withdraw.noae <- beta.par(0.405000000)$beta
  
  # The standard deviation for the data without adverse event
  sig.X.noae <- 0.300
  
  # The standard deviation for withdrawn
  sig.X.with <- 0.310

  # Random generation for the Beta distribution with parameters for utility of 
    # no adverse effect
  X2 <- rbeta(n,
              betaPar(theta$u.l1.noae, sig.X.noae)$a,
              betaPar(theta$u.l1.noae, sig.X.noae)$b)
  
  # Random generation for the Beta distribution with parameters for utility 
    # of adverse effect
  X4 <- rbeta(n,
              betaPar(theta$u.l1.withdraw.noae, sig.X.with)$a,
              betaPar(theta$u.l1.withdraw.noae, sig.X.with)$b)
  
  # Random generation for binomial distribution with parameters size and probability
  missingness <- rbinom(n, 1, 0.687)
  
  X2[which(X2*missingness==0)] <- NA
  X4[which(X4*missingness==0)] <- NA
  
  # Summarize simulated data as geometric mean
  X.gam <- c(psych::geometric.mean(X2),
             psych::geometric.mean(1-X2))
  X.gam.with <- c(psych::geometric.mean(X4),
                  psych::geometric.mean(1-X4))
  
  # Output simulated data and its summary statistics
  r <- list()
  
  # Sum(rand) x 2 matrices for Beta distribution parameters 
  r$simulated_data <- list(r.u.l1.noae, 
                           r.u.l1.withdraw.noae,
                           s.u.l1.noae, 
                           s.u.l1.withdraw.noae,
                           X2, 
                           X4,
                           sig.X.noae, 
                           sig.X.with)
  
  names(r$simulated_data) <- c("r.u.l1.noae", 
                               "r.u.l1.withdraw.noae",
                               "s.u.l1.noae", 
                               "s.u.l1.withdraw.noae",
                               "utility.noae", 
                               "utility.ae",
                               "sig.X.noae", 
                               "sig.X.with")
  
  # Sum(rand) x 2 matrices for summarized statistics
  r$summary_statistics <- list("X.gam.1"=X.gam[1], 
                               "X.gam.2"=X.gam[2],
                               "X.gam.with.1"=X.gam.with[1], 
                               "X.gam.with.2"=X.gam.with[2])
  
  return(r)
}

likelihood_compute_function <- function(n, D, theta){
  l <- array(NA, dim=nrow(theta))
  
  sig.X.noae=0.300
  sig.X.with=0.310
  
  # Beta distribution parameter for the utility of no adverse effects
  noae.par <- betaPar(theta$u.l1.noae, sig.X.noae)
  
  # Beta distribution parameter for the utility of adverse effects
  ae.par <- betaPar(theta$u.l1.ae, sig.X.with)
  
  # Likelihood function for the data
  for(k in 1:nrow(theta)) {
    l[k] <- exp(sum(dbeta(D[,1], noae.par$a[k], noae.par$b[k], log = T), na.rm=T) +
                  sum(dbeta(D[,2], ae.par$a[k], ae.par$b[k], log = T), na.rm=T))
  }
  
  # Set NA and Infinite values to zero
  l[is.na(l)] <- 0
  l[is.infinite(l)] <- 0
  
  return(l)
}

evsi_specification <- list(evppi_parameters=c("u.l1.noae","u.l1.withdraw.noae"),
                           gam_tensor_prod_dim=3,
                           likelihood_data_variables=c("utility.noae", "utility.ae"),
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

decision_model_cea <- list(cost=discount.15*cbind(costs.t1,costs.t2),
                           effect=discount.15*cbind(utility.t1,utility.t2),
                           wtp=20000)

decision_model_psa <- as.data.frame(pars)

data_generation_function <- function(n, theta){
  #Utilities
  s.u.l1.noae <- beta.par(0.695000000)$alpha
  r.u.l1.noae <- beta.par(0.695000000)$beta
  
  #Withdraw due to other reasons
  s.u.l1.withdraw.noae <- beta.par(0.405000000)$alpha
  r.u.l1.withdraw.noae <- beta.par(0.405000000)$beta
  
  # The standard deviation for the data without adverse event
  sig.X.noae <- 0.300
  
  # The standard deviation for withdrawn
  sig.X.with <- 0.310
  
  # Random generation for the Beta distribution with parameters for utility 
    # of no adverse effect
  X2 <- rbeta(n,
              betaPar(theta$u.l1.noae, sig.X.noae)$a,
              betaPar(theta$u.l1.noae, sig.X.noae)$b)
  
  # Random generation for the Beta distribution with parameters for utility 
    # of adverse effect
  X4 <- rbeta(n,
              betaPar(theta$u.l1.withdraw.noae, sig.X.with)$a,
              betaPar(theta$u.l1.withdraw.noae, sig.X.with)$b)
  
  # Random generation for binomial distribution with parameters size and probability
  missingness <- rbinom(n, 1, 0.687)
  
  X2[which(X2*missingness==0)] <- NA
  X4[which(X4*missingness==0)] <- NA
  
  # Summarize simulated data as geometric mean
  X.gam <- c(psych::geometric.mean(X2),
             psych::geometric.mean(1-X2))
  X.gam.with <- c(psych::geometric.mean(X4),
                  psych::geometric.mean(1-X4))
  
  # Output simulated data and its summary statistics
  r <- list()
  
  # Sum(rand) x 2 matrices for Beta distribution parameters 
  r$simulated_data <- list(r.u.l1.noae, 
                           r.u.l1.withdraw.noae,
                           s.u.l1.noae, 
                           s.u.l1.withdraw.noae,
                           X2, 
                           X4,
                           sig.X.noae, 
                           sig.X.with)
  
  names(r$simulated_data) <- c("r.u.l1.noae", 
                               "r.u.l1.withdraw.noae",
                               "s.u.l1.noae", 
                               "s.u.l1.withdraw.noae",
                               "utility.noae", 
                               "utility.ae",
                               "sig.X.noae", 
                               "sig.X.with")
  
  # Sum(rand) x 2 matrices for summarized statistics
  r$summary_statistics <- list("X.gam.1"=X.gam[1], 
                               "X.gam.2"=X.gam[2],
                               "X.gam.with.1"=X.gam.with[1], 
                               "X.gam.with.2"=X.gam.with[2])
  
  return(r)
}

evsi_specification <- list(gam_parameters=c("X.gam.1", 
                                            "X.gam.2",
                                            "X.gam.with.1", 
                                            "X.gam.with.2"),
                           gam_tensor_prod_dim=3,
                           outer_cond_exp_loop_size=10000,
                           sim_data_size=c(150))

evsi_regression_based(decision_model_cea,
                      decision_model_psa,
                      data_generation_function,
                      evsi_specification)
