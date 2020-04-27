###################################################################################
# Ades et al
###################################################################################

############################################################
# Load model and EVSI methods scripts
############################################################
source(".//Decision Model Example - Ades et al//Ades_et_al_decision_model.R")
source(".//EVSI Estimation Methods//evsi_moment_matching_method.R")
source(".//EVSI Estimation Methods//evsi_Gaussian_approximation_method.R")
source(".//EVSI Estimation Methods//evsi_importance_sampling_method.R")
source(".//EVSI Estimation Methods//evsi_regression_based_method.R")

############################################################
# Moment Matching Method
############################################################

cost_funtion <- function(L, Qe, Qse, Ce, Ct, Cse, Pc, Pse, Pt, lambda) {
  cost.standard <- 
    Pc*(L*(1+Qe)/2)+
    (1-Pc)*L
  
  cost.treatment <- 
    Pse*Pt*((L*(1+Qe)/2-Qse))+
    Pse*(1-Pt)*((L-Qse))+
    (1-Pse)*Pt*(L*(1+Qe)/2)+
    (1-Pse)*(1-Pt)*(L)
  
  return(c(cost.standard, cost.treatment))
}

effect_funtion <- function(L, Qe, Qse, Ce, Ct, Cse, Pc, Pse, Pt, lambda) {
  utility.standard <- 
    Pc*(Ce)+
    (1-Pc)*0
  
  utility.treatment <- 
    Pse*Pt*((Ct+Cse+Ce))+
    Pse*(1-Pt)*((Ct+Cse))+
    (1-Pse)*Pt*((Ct+Ce))+
    (1-Pse)*(1-Pt)*(Ct)
  
  return(c(utility.standard, utility.treatment))
}

decision_model_cea <- list(cost=as.matrix(cost),
                           effect=as.matrix(utility),
                           cost_funtion=cost_funtion,
                           effect_funtion=effect_funtion,
                           wtp=75000)

decision_model_psa <- theta

# Data Collection Scenario 1: EVSI for the Probability of Side Effects (Pse)

model.data.full <- function(){
  L<-30
  Qse<-1
  Ce<-200000
  Ct<-15000
  Cse<-100000
  Pc~dbeta(15,85)
  Pse~dbeta(3,9)
  OR.log~dnorm(-1.5,3)
  OR<-exp(OR.log)
  Pt<-exp(log(Pc/(1-Pc))+log(OR))/(exp(log(Pc/(1-Pc))+log(OR))+1)
  lambda<-75000
  
  X~dbin(Pse,60)
  lQe~dnorm(0.6,6)
  #X~dnorm(lQe,1000000)
  Qe<-exp(lQe)/(exp(lQe)+1)
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=list(),
                            inits=NULL,
                            iter_number=10000,
                            burnin_number=1000,
                            mcmc_parameter=c("Pse"))

data_generation_function <- function(n, theta){
  X <- rbinom(1, n, theta$Pse)
  Pse.post <- mean(rbeta(1000, 3+X, 69-X))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("X"=X, "Pse"=Pse.post)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(evppi_parameters=c("Pse"),
                           nested_sample_size=50,
                           quantile_parameter=c("Pse"),
                           sim_data_size=c(60))

evsi_moment_matching(decision_model_cea,
                     decision_model_psa,
                     decision_model_mcmc,
                     data_generation_function,
                     evsi_specification)

# Data Collection Scenario 2: EVSI for Quality of Life after Critical Event (Qe)

model.data.full <- function(){
  L<-30
  Qse<-1
  Ce<-200000
  Ct<-15000
  Cse<-100000
  Pc~dbeta(15,85)
  Pse~dbeta(3,9)
  OR.log~dnorm(-1.5,3)
  OR<-exp(OR.log)
  Pt<-exp(log(Pc/(1-Pc))+log(OR))/(exp(log(Pc/(1-Pc))+log(OR))+1)
  lambda<-75000
  
  #X~dbin(Pse,60)
  lQe~dnorm(0.6,6)
  X~dnorm(lQe,50)
  Qe<-exp(lQe)/(exp(lQe)+1)
  
  NB1<-Pc*(lambda*L*(1+Qe)/2-Ce)+
    (1-Pc)*lambda*L
  
  NB2<-Pse*Pt*(lambda*(L*(1+Qe)/2-Qse)-(Ct+Cse+Ce))+
    Pse*(1-Pt)*(lambda*(L-Qse)-(Ct+Cse))+
    (1-Pse)*Pt*(lambda*L*(1+Qe)/2-(Ct+Ce))+
    (1-Pse)*(1-Pt)*(lambda*L-Ct)
  INB<-NB2-NB1
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=list(),
                            update_prior_parameter=c("mu.X","lQe"),
                            inits=NULL,
                            iter_number=10000,
                            burnin_number=1000,
                            mcmc_parameter=c("Qe"))

data_generation_function <- function(n, theta){
  X <- rnorm(n, mean=logit(theta$Qe), sd=sqrt(1/50))
  lQe.p <- rnorm(1000, (0.6*6+sum(X)*50)/(n*50+6), sqrt(1/(n*50+6)))
  lQe.post <- mean(lQe.p)
  Qe.post <- mean(inv.logit(lQe.p))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("X"=X, 
                           "Qe"=Qe.post, 
                           "lQe"=lQe.post,
                           "mu.X"=mean(X))
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(evppi_parameters=c("Qe"),
                           nested_sample_size=50,
                           quantile_parameter=c("Qe"),
                           sim_data_size=c(100))

evsi_moment_matching(decision_model_cea,
                     decision_model_psa,
                     decision_model_mcmc,
                     data_generation_function,
                     evsi_specification)

# Data Collection Scenario 3: EVSI for the Treatment Effect Size (OR)

model.data.full <- function(){
  L<-30
  Qse<-1
  Ce<-200000
  Ct<-15000
  Cse<-100000
  Pc.post~ dbeta(15,85)
  Pc~dbeta(15,85)
  Pse~dbeta(3,9)
  OR.log~dnorm(-1.5,3)
  OR<-exp(OR.log)
  logit(Pt.post)=logit(Pc.post)+OR.log
  logit(Pt)<-logit(Pc)+OR.log
  #Pt<-exp(log(Pc/(1-Pc))+log(OR))/(exp(log(Pc/(1-Pc))+log(OR))+1)
  lambda<-75000
  
  Dc~dbin(Pc.post,200) 
  Dt~dbin(Pt.post,200)
  
  lQe~dnorm(0.6,6)
  Qe<-exp(lQe)/(exp(lQe)+1)
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=list(),
                            update_prior_parameter=c("Dc","Dt"),
                            inits=NULL,
                            iter_number=10000,
                            burnin_number=1000,
                            mcmc_parameter=c("OR"))

data_generation_function <- function(n, theta){
  Dc <- rbinom(1, n, theta$Pc)
  Dt <- rbinom(1, n, theta$Pt)
  
  var.OR <- 1/Dc+1/(200-Dc)+1/Dt+1/(200-Dt)
  mean.OR <- log(Dt/(200-Dt)/(Dc/(200-Dc)))
  if(Dc==0 | Dt==0){
    var.OR <- 1.5/(1+Dc)+1/(200-Dc)+1.5/(1+Dt)+1/(200-Dt)
    mean.OR <- log((Dt+0.5)/(201-Dt)/((Dc+0.5)/(201-Dc)))
  }
  mean.pr <- -1.5
  var.pr <- 1/3
  var.post <- 1/(1/var.OR+3)
  mean.post <- (1/(var.OR)*mean.OR+mean.pr*3)*var.post
  OR.post <- mean(exp(rnorm(1000, mean.post, sqrt(var.post))))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("Dc"=Dc, "Dt"=Dt,
                           "OR"=OR.post)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(evppi_parameters=c("OR"),
                           nested_sample_size=50,
                           quantile_parameter=c("Pc", "Pt"),
                           sim_data_size=c(100))

evsi_moment_matching(decision_model_cea,
                     decision_model_psa,
                     decision_model_mcmc,
                     data_generation_function,
                     evsi_specification)

############################################################
# Gaussian Approximation Method
############################################################

decision_model_cea <- list(cost=as.matrix(cost),
                           effect=as.matrix(utility),
                           wtp=75000)

decision_model_psa <- theta

# Data Collection Scenario 1: EVSI for the Probability of Side Effects (Pse)
model.data.full <- function(){
  L<-30
  Qse<-1
  Ce<-200000
  Ct<-15000
  Cse<-100000
  Pc~dbeta(15,85)
  Pse~dbeta(3,9)
  OR.log~dnorm(-1.5,3)
  OR<-exp(OR.log)
  Pt<-exp(log(Pc/(1-Pc))+log(OR))/(exp(log(Pc/(1-Pc))+log(OR))+1)
  lambda<-75000
  
  X~dbin(Pse,60)
  lQe~dnorm(0.6,6)
  #X~dnorm(lQe,1000000)
  Qe<-exp(lQe)/(exp(lQe)+1)
  NB1<-Pc*(lambda*L*(1+Qe)/2-Ce)+
    (1-Pc)*lambda*L
  
  NB2<-Pse*Pt*(lambda*(L*(1+Qe)/2-Qse)-(Ct+Cse+Ce))+
    Pse*(1-Pt)*(lambda*(L-Qse)-(Ct+Cse))+
    (1-Pse)*Pt*(lambda*L*(1+Qe)/2-(Ct+Ce))+
    (1-Pse)*(1-Pt)*(lambda*L-Ct)
  INB<-NB2-NB1
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=list(),
                            update_prior_parameter=c("X","Pse"),
                            inits=NULL,
                            chain_number=3,
                            burnin_number=1000,
                            mcmc_parameter=c("Pse"))

data_generation_function <- function(n, theta){
  X <- rbinom(1, n, theta$Pse)
  Pse.post <- mean(rbeta(1000, 3+X, 69-X))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("X"=X, "Pse"=Pse.post)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(evppi_parameters=c("Pse"),
                           economic_model_size=60,
                           outer_cond_exp_loop_size=10000,
                           inner_cond_exp_loop_size=5000,
                           sim_data_size=c(60),
                           mcmc_retry=5)

evsi_gaussian_approximation(decision_model_cea,
                            decision_model_psa,
                            decision_model_mcmc,
                            data_generation_function,
                            evsi_specification)

# Data Collection Scenario 2: EVSI for Quality of Life after Critical Event (Qe)
model.data.full <- function(){
  L<-30
  Qse<-1
  Ce<-200000
  Ct<-15000
  Cse<-100000
  Pc~dbeta(15,85)
  Pse~dbeta(3,9)
  OR.log~dnorm(-1.5,3)
  OR<-exp(OR.log)
  Pt<-exp(log(Pc/(1-Pc))+log(OR))/(exp(log(Pc/(1-Pc))+log(OR))+1)
  lambda<-75000
  
  #X~dbin(Pse,60)
  lQe~dnorm(0.6,6)
  X~dnorm(lQe,50)
  Qe<-exp(lQe)/(exp(lQe)+1)
  
  NB1<-Pc*(lambda*L*(1+Qe)/2-Ce)+
    (1-Pc)*lambda*L
  
  NB2<-Pse*Pt*(lambda*(L*(1+Qe)/2-Qse)-(Ct+Cse+Ce))+
    Pse*(1-Pt)*(lambda*(L-Qse)-(Ct+Cse))+
    (1-Pse)*Pt*(lambda*L*(1+Qe)/2-(Ct+Ce))+
    (1-Pse)*(1-Pt)*(lambda*L-Ct)
  INB<-NB2-NB1
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=list(),
                            update_prior_parameter=c("mu.X","lQe"),
                            inits=NULL,
                            chain_number=3,
                            burnin_number=1000,
                            mcmc_parameter=c("Qe"))

data_generation_function <- function(n, theta){
  X <- rnorm(n, mean=logit(theta$Qe), sd=sqrt(1/50))
  
  lQe.p <- rnorm(1000, (0.6*6+sum(X)*50)/(n*50+6), sqrt(1/(n*50+6)))
  lQe.post <- mean(lQe.p)
  Qe.post <- mean(inv.logit(lQe.p))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("X"=X, "Qe"=Qe.post, "lQe"=lQe.post,
                           "mu.X"=mean(X))
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(evppi_parameters=c("Qe"),
                           economic_model_size=100,
                           outer_cond_exp_loop_size=10000,
                           inner_cond_exp_loop_size=5000,
                           sim_data_size=c(100),
                           mcmc_retry=5)

evsi_gaussian_approximation(decision_model_cea,
                            decision_model_psa,
                            decision_model_mcmc,
                            data_generation_function,
                            evsi_specification)

# Data Collection Scenario 3: EVSI for the Treatment Effect Size (OR)
model.data.full <- function(){
  L<-30
  Qse<-1
  Ce<-200000
  Ct<-15000
  Cse<-100000
  Pc.post~ dbeta(15,85)
  Pc~dbeta(15,85)
  Pse~dbeta(3,9)
  OR.log~dnorm(-1.5,3)
  OR<-exp(OR.log)
  logit(Pt.post)=logit(Pc.post)+OR.log
  logit(Pt)<-logit(Pc)+OR.log
  #Pt<-exp(log(Pc/(1-Pc))+log(OR))/(exp(log(Pc/(1-Pc))+log(OR))+1)
  lambda<-75000
  
  Dc~dbin(Pc.post,200) 
  Dt~dbin(Pt.post,200)
  
  lQe~dnorm(0.6,6)
  Qe<-exp(lQe)/(exp(lQe)+1)
  
  NB1<-Pc*(lambda*L*(1+Qe)/2-Ce)+
    (1-Pc)*lambda*L
  
  NB2<-Pse*Pt*(lambda*(L*(1+Qe)/2-Qse)-(Ct+Cse+Ce))+
    Pse*(1-Pt)*(lambda*(L-Qse)-(Ct+Cse))+
    (1-Pse)*Pt*(lambda*L*(1+Qe)/2-(Ct+Ce))+
    (1-Pse)*(1-Pt)*(lambda*L-Ct)
  INB<-NB2-NB1
}

decision_model_mcmc <- list(model=model.data.full,
                            prior_data=list(),
                            update_prior_parameter=c("Dc","Dt"),
                            inits=NULL,
                            chain_number=3,
                            burnin_number=1000,
                            mcmc_parameter=c("OR"))

data_generation_function <- function(n, theta){
  Dc <- rbinom(1, n, theta$Pc)
  Dt <- rbinom(1, n, theta$Pt)
  
  var.OR <- 1/Dc+1/(200-Dc)+1/Dt+1/(200-Dt)
  mean.OR <- log(Dt/(200-Dt)/(Dc/(200-Dc)))
  if(Dc==0 | Dt==0){
    var.OR <- 1.5/(1+Dc)+1/(200-Dc)+1.5/(1+Dt)+1/(200-Dt)
    mean.OR <- log((Dt+0.5)/(201-Dt)/((Dc+0.5)/(201-Dc)))
  }
  mean.pr <- -1.5
  var.pr <- 1/3
  var.post <- 1/(1/var.OR+3)
  mean.post <- (1/(var.OR)*mean.OR+mean.pr*3)*var.post
  OR.post <- mean(exp(rnorm(1000, mean.post, sqrt(var.post))))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("Dc"=Dc, "Dt"=Dt,
                           "OR"=OR.post)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(evppi_parameters=c("OR"),
                           economic_model_size=200,
                           outer_cond_exp_loop_size=10000,
                           inner_cond_exp_loop_size=5000,
                           sim_data_size=c(200),
                           mcmc_retry=5)

evsi_gaussian_approximation(decision_model_cea,
                            decision_model_psa,
                            decision_model_mcmc,
                            data_generation_function,
                            evsi_specification)

############################################################
# Importance Sampling Method
############################################################

decision_model_cea <- list(cost=cost,
                           effect=utility,
                           wtp=75000)

decision_model_psa <- as.data.frame(theta)

# Data Collection Scenario 1: EVSI for the Probability of Side Effects (Pse)
data_generation_function <- function(n, theta){
  X <- rbinom(1, n, theta$Pse)
  Pse.post <- mean(rbeta(1000, 3+X, 69-X))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("X"=X, "Pse"=Pse.post)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

likelihood_compute_function <- function(n, D, theta){

  l <- dbinom(D[,"X"], n, theta[,"Pse"], log = TRUE)
  
  return(exp(l))
}

evsi_specification <- list(evppi_parameters=c("Pse"),
                           gam_tensor_prod_dim=NA,
                           likelihood_data_variables=c("X"),
                           outer_cond_exp_loop_size=10000,
                           sim_data_size=c(60))

evsi_importance_sample(decision_model_cea,
                       decision_model_psa,
                       data_generation_function,
                       likelihood_compute_function,
                       evsi_specification)

# Data Collection Scenario 2: EVSI for Quality of Life after Critical Event (Qe)

data_generation_function <- function(n, theta){
  X <- rnorm(n, mean=logit(theta$Qe), sd=sqrt(1/50))
  
  lQe.p <- rnorm(1000, (0.6*6+sum(X)*50)/(n*50+6), sqrt(1/(n*50+6)))
  lQe.post <- mean(lQe.p)
  Qe.post <- mean(inv.logit(lQe.p))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("X"=X, "Qe"=Qe.post, "lQe"=lQe.post,
                           "mu.X"=mean(X))
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

likelihood_compute_function <- function(n, D, theta){
  l <- array(NA, dim=nrow(theta))
  
  for(k in 1:nrow(theta)) {
    l[k] <- sum(dnorm(D[,"X"], mean=logit(theta[k,"Qe"]), sd=sqrt(1/50), 
                      log = TRUE))
  }
  
  return(exp(l))
}

evsi_specification <- list(evppi_parameters=c("Qe"),
                           gam_tensor_prod_dim=NA,
                           likelihood_data_variables=c("X"),
                           outer_cond_exp_loop_size=10000,
                           sim_data_size=c(100))

evsi_importance_sample(decision_model_cea,
                       decision_model_psa,
                       data_generation_function,
                       likelihood_compute_function,
                       evsi_specification)

# Data Collection Scenario 3: EVSI for the Treatment Effect Size (OR)

data_generation_function <- function(n, theta){
  Dc <- rbinom(1, n, theta$Pc)
  Dt <- rbinom(1, n, theta$Pt)

  var.OR <- 1/Dc+1/(200-Dc)+1/Dt+1/(200-Dt)
  mean.OR <- log(Dt/(200-Dt)/(Dc/(200-Dc)))
  if(Dc==0 | Dt==0){
    var.OR <- 1.5/(1+Dc)+1/(200-Dc)+1.5/(1+Dt)+1/(200-Dt)
    mean.OR <- log((Dt+0.5)/(201-Dt)/((Dc+0.5)/(201-Dc)))
  }
  mean.pr <- -1.5
  var.pr <- 1/3
  var.post <- 1/(1/var.OR+3)
  mean.post <- (1/(var.OR)*mean.OR+mean.pr*3)*var.post
  OR.post <- mean(exp(rnorm(1000, mean.post, sqrt(var.post))))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("Dc"=Dc, "Dt"=Dt,
                           "OR"=OR.post)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

likelihood_compute_function <- function(n, D, theta){
  
  l <- dbinom(D[,"Dc"], n, theta[,"Pc"], log = TRUE) +
    dbinom(D[,"Dt"], n, theta[,"Pt"], log = TRUE)

  return(exp(l))
}

evsi_specification <- list(evppi_parameters=c("OR"),
                           gam_tensor_prod_dim=NA,
                           likelihood_data_variables=c("Dc", "Dt"),
                           outer_cond_exp_loop_size=10000,
                           sim_data_size=c(200))

evsi_importance_sample(decision_model_cea,
                       decision_model_psa,
                       data_generation_function,
                       likelihood_compute_function,
                       evsi_specification)

############################################################
# Regression-Based Method
############################################################

decision_model_cea <- list(cost=cost,
                           effect=utility,
                           wtp=75000)

decision_model_psa <- theta

# Data Collection Scenario 1: EVSI for the Probability of Side Effects (Pse)

data_generation_function <- function(n, theta){
  X <- rbinom(1, n, theta$Pse)
  Pse.post <- mean(rbeta(1000, 3+X, 69-X))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("X"=X, "Pse"=Pse.post)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(gam_parameters=c("Pse"),
                           gam_tensor_prod_dim=NA,
                           outer_cond_exp_loop_size=10000,
                           sim_data_size=c(60))

evsi_regression_based(decision_model_cea,
                      decision_model_psa,
                      data_generation_function,
                      evsi_specification)

# Data Collection Scenario 2: EVSI for Quality of Life after Critical Event (Qe)

data_generation_function <- function(n, theta){
  X <- rnorm(n, mean=logit(theta$Qe), sd=sqrt(1/50))
  
  lQe.p <- rnorm(1000, (0.6*6+sum(X)*50)/(n*50+6), sqrt(1/(n*50+6)))
  lQe.post <- mean(lQe.p)
  Qe.post <- mean(inv.logit(lQe.p))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("X"=X, "Qe"=Qe.post, "lQe"=lQe.post,
                           "mu.X"=mean(X))
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(gam_parameters=c("Qe"),
                           gam_tensor_prod_dim=NA,
                           outer_cond_exp_loop_size=10000,
                           sim_data_size=c(100))

evsi_regression_based(decision_model_cea,
                      decision_model_psa,
                      data_generation_function,
                      evsi_specification)

# Data Collection Scenario 3: EVSI for the Treatment Effect Size (OR)

data_generation_function <- function(n, theta){
  Dc <- rbinom(1, n, theta$Pc)
  Dt <- rbinom(1, n, theta$Pt)
  
  var.OR <- 1/Dc+1/(200-Dc)+1/Dt+1/(200-Dt)
  mean.OR <- log(Dt/(200-Dt)/(Dc/(200-Dc)))
  if(Dc==0 | Dt==0){
    var.OR <- 1.5/(1+Dc)+1/(200-Dc)+1.5/(1+Dt)+1/(200-Dt)
    mean.OR <- log((Dt+0.5)/(201-Dt)/((Dc+0.5)/(201-Dc)))
  }
  mean.pr <- -1.5
  var.pr <- 1/3
  var.post <- 1/(1/var.OR+3)
  mean.post <- (1/(var.OR)*mean.OR+mean.pr*3)*var.post
  OR.post <- mean(exp(rnorm(1000, mean.post, sqrt(var.post))))
  
  # Output simulated data and its summary statistics
  D <- list()
  
  D$simulated_data <- list("Dc"=Dc, "Dt"=Dt,
                           "OR"=OR.post)
  
  # No summary statistic transformation applied
  D$summary_statistics <- D$simulated_data
  
  return(D)
}

evsi_specification <- list(gam_parameters=c("OR"),
                           gam_tensor_prod_dim=NA,
                           outer_cond_exp_loop_size=10000,
                           sim_data_size=c(200))

evsi_regression_based(decision_model_cea,
                      decision_model_psa,
                      data_generation_function,
                      evsi_specification)
