###################################################################################
#  # A brief step-by-step guide for Gaussian Approximation Method #
#  1. PSA simulation to obtain theta_s and INB_t(theta_s), s = 1,...,S. 
#     (Performed by running chemotherapy_model.R).
#  2. Fit T-1 regression models with INB_t(theta_s) as outcomes and phi_s as covariates.
#  3. For each of the P elements in phi:
#     (a) Determine the prior effective sample size N_0(p). Methods to estimate 
#         the prior effective sample size are presented in the supplementary material.
#     (b) For a proposed study with N participants, compute a weighted sum of 
#         phi and phi_bar (the mean of the parameters) by rescaling the simulations 
#         for the pth element of phi by multiplying by N/(N+N_0(p)) and 
#         multiplying phi_bar by 1-[N/(N+N_0(p))].
#  4. Using the regression models from step 2, predict the model outcomes for
#     the rescaled phi simulations.
#  5. The fitted values from step 4 estimate mu_t(X_s), which is the EVSI. 
###################################################################################

############################################################
# Function to run the Gaussian Approximation Method
############################################################
evsi_gaussian_approximation <- functionfunction(
                                 # Cost effective analysis for the decision model
                                 decision_model_cea,
                                 # Matrix of parameters of interest from baseline model
                                 decision_model_psa,
                                 # Information relevant to MCMC simulation to resample parameters
                                 decision_model_mcmc,
                                 # Simulate the dataset X_s from p(X|phi_s)
                                 data_generation_function,
                                 # Parameters of interest for the model
                                 evsi_specification) {
  
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
  
  p_install(earth, force=FALSE)       # Multivariate Adaptive Regression Splines
  p_load(earth)
  
  p_install(R2jags, force=FALSE)      # Using R to run 'JAGS'
  p_load(R2jags)
  
  p_install(R2OpenBUGS, force=FALSE)  # Using R to run 'OpenBUGS'
  p_load(R2OpenBUGS)
  
  p_install(BCEA, force=FALSE)        # Bayesian cost-effectiveness Analysis
  p_load(BCEA)
  
  # Mixed GAM computation vehicle with automatic smoothness estimation
  p_install(mgcv, force=FALSE)        
  p_load(mgcv)
  
  # Expected value of sample information
  if(!require(EVSI)){  
    # To install, use devtools::install_github("annaheath/EVSI")
    devtools::install_github("annaheath/EVSI")  
    library(EVSI)
  }
  
  ############################################################
  # Load required helper script and funciton
  ############################################################
  source(choose.files(caption="Select script containing Gaussian approximation function (default: predict_ga.R)",
                      multi=F), encoding='WINDOWS-1252')
  
  ############################################################
  # Baseline cost and effectiveness analysis
  ############################################################
  
  # Performs baseline cost-effectiveness analysis for specified willingness to pay
  NB <- decision_model_cea$effect * decision_model_cea$wtp - decision_model_cea$cost
  
  # Find incremental net benefit for each treatment
  d.star <- which.max(apply(NB, 2, mean))
  INB <- NB - NB[, d.star]
  INB.no_reference <- INB[, -d.star, drop=F]

  ############################################################
  # Estimate EVSI by the Gaussian Approximation Method
  ############################################################
  
  ### Step1: PSA simulation to obtain theta_s and INB_t(theta_s)
  
  ### Step2: Fit T-1 regression model with INB(y) as outcomes and theta as covariates

  # Fit T-1 regression model
  gam_formula <- paste0("INB_t ~ s(",
                        paste0(evsi_specification$evppi_parameter, collapse=")+s("),
                        ")")
  
  gam_model <- list()
  
  for(k in 1:ncol(INB.no_reference)) {
    INB_t <- INB.no_reference[,k]
    gam_model[[k]] <- gam(as.formula(gam_formula), data=decision_model_psa)
  }

  ### Step3(a): Determine the prior effective sample size N_0(p)
  
  # Specify an economic model with n
  n <- evsi_specification$economic_model_size
  
  # Set the number of chains in iterations
  Size.Outer <- evsi_specification$outer_cond_exp_loop_size
  Size.Inner <- evsi_specification$inner_cond_exp_loop_size
  
  # Number of chains
  n.chains <- decision_model_mcmc$chain_number  
  # Number of burn in iterations
  n.burnin <- decision_model_mcmc$burnin_number    
  # Number of iterations per chain
  n.iter <- ceiling(Size.Inner/n.chains) + n.burnin    
  
  # Generate Bayesian Model File 
  filein <- file.path(tempdir(), fileext="datmodel.txt")
  write.model(decision_model_mcmc$model, filein)
  
  # Generate data
  for(i in 1:Size.Outer){
    mcmc_pass <- FALSE
    mcmc_run <- 1
    
    while(!mcmc_pass & (mcmc_run <= evsi_specification$mcmc_retry)) {
      # Simulate dataset
      sim_data <- data_generation_function(n=evsi_specification$economic_model_size,
                                           theta=decision_model_psa[i,])
      
      # Prior data and parameters for JAGS model
      data.full <- append(decision_model_mcmc$prior_data, 
                          sim_data$simulated_data[decision_model_mcmc$update_prior_parameter])
      data.full <- append(data.full, list("n"=n))
      
      # Perform the MCMC simulation with Jags
      tryCatch({
        jags.data <- jags(data=data.full,
                          inits=decision_model_mcmc$inits,
                          parameters.to.save=decision_model_mcmc$mcmc_parameter,
                          model.file=filein,
                          n.chains=n.chains,
                          n.iter=n.iter,
                          n.thin=1,
                          n.burnin=n.burnin)
        
        if(exists("mean.data")) {
          mean.data <- rbind(mean.data,
                             apply(as.data.frame(jags.data$BUGSoutput$sims.list),
                                   c(2), mean))
        } else {
          mean.data <- apply(as.data.frame(jags.data$BUGSoutput$sims.list),
                             c(2), mean)
        }
        
        mcmc_pass <- TRUE
        
      }, error=function(e, mcmc_run.=mcmc_run) {
        print(paste0("Invalid MCMC: ", mcmc_run.))
      })
      
      # Number of MCMC run for each PSA sample
      mcmc_run <- mcmc_run + 1
    }
    
  }
  
  # Estimating n0 separately for each element of the meta-model
  # Determine the scalar or vector of effective prior size n0 
  var.theta <- apply(decision_model_psa[,colnames(decision_model_psa) %in% colnames(mean.data),drop=F],
                     c(2), var)

  n0 <- lapply(names(var.theta), function(x, n, var.theta, mean.data) {
    return(n * (var.theta[x] / var(mean.data[,x]) - 1))
  }, n=n, var.theta=var.theta, mean.data=mean.data)
  
  n0 <- unlist(n0)
  
  ### Step3(b): For a proposed study with N participants, compute a weighted sum 
                # of phi and phi_bar (the mean of the parameters) by rescaling 
                # the simulations for the pth element of phi by multiplying by 
                # N/(N+N_0(p)) and multiplying phi_bar by 1-[N/(N+N_0(p))]
  
  evsi.jalal <- array(NA, dim=length(evsi_specification$sim_data_size))
  
  for(i in 1:length(evsi_specification$sim_data_size)) {
    print(paste0("EVSI estimation ", i, 
                 ": simulated data n=", evsi_specification$sim_data_size[i]))
    
    ### Step4: Using the regression models from step 2, predict the model 
               # outcomes for the rescaled phi simulations
    llpred.n <- rep(evsi_specification$sim_data_size[i],
                    # scalar or vector of new sample size to compute evsi on
                    length(n0))  

    llpred <- sapply(1:length(gam_model), function(k) {
      return(predict.ga(gam_model[[k]], n0=n0, n=llpred.n))
    })
    
    ### Step5: The fitted values from step 4 estimate mu_t(X_s), which is the EVSI
    llpred_max <- apply(llpred, c(1), max)
    
    # Estimates EVSI
    evsi.jalal[i] <- mean(pmax(0, llpred_max)) - max(0, mean(llpred_max))
  }
  
  return(evsi.jalal)
}


 
