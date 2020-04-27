###################################################################################
#  # A brief step-by-step guide for the Moment Matching Method #
#  1. PSA simulation to obtain theta_s and INB_t(theta_s) 
#     (Performed by running chemotherapy_model.R);
#  2. Estimate EVPPI to obtain simulations of the inner conditional expectation
#     of INB on the future data (eta_t(s));
#  3. Extract Q, which is the outer loop size as well as the number of nested 
#     samples (30 < Q < 50), sample quantiles from the simulations of the 
#     parameters of interest (phi);
#  4. For q = 1,...,Q:
#     (a) Simulate a future dataset from the sampling distribution p(X|phi)
#     (b) Use Bayesian methods to update the distribution of the model parameters
#     (c) Rerun the probabilistic sensitivity analysis to update the distribution
#         of the net monetary benefit
#     (d) Calculate the variance of the net monetary benefit "sigma_q"
#  5. Calculate the sigma^2 = Var[INB_t(theta_s)] - (1/Q)*sum(sigma_q^2)
#  6. Rescale the simulations of eta_t(s) such that their variance is equal to sigma^2;
#  7. These rescaled simulations estimate mu_t(X_s), which is the EVSI 
###################################################################################

############################################################
# Function to run the Moment Matching Method
############################################################
evsi_moment_matching <- function(# Cost effective analysis for the decision model
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
  # Baseline cost and effectiveness analysis
  ############################################################
  bcea_result <- bcea(decision_model_cea$effect,
                      decision_model_cea$cost,
                      ref=2, wtp=decision_model_cea$wtp)
  
  ############################################################
  # Estimate EVSI by Moment Matching Method
  ############################################################
  
  ### Step1: PSA simulation to obtain theta_s and INB_t(theta_s) 
             # (Performed by running chemotherapy_model.R)
 
  ### Step2: Estimate EVPPI to obtain simulations of the inner conditional 
             # expectation of INB on the future data (eta_t(s))
  evppi.full <- evppi(evsi_specification$evppi_parameter,
                      decision_model_psa,
                      bcea_result)
  
  ### Step3: Extract Q which is the number of nested samples and sample quantiles 
             # from the simulations of the parameters of interest (phi)
  Q <- evsi_specification$nested_sample_size
  
  # Estimate EVSI for different sample size
  evsi.heath <- array(NA, dim=length(evsi_specification$sim_data_size))
  
  for(i in 1:length(evsi_specification$sim_data_size)) {
    print(paste0("EVSI estimation ", i, 
                 ": simulated data n=", evsi_specification$sim_data_size[i]))
    
    # Generate the future data from specified sampling distributions
    quants <- gen.quantiles(evsi_specification$quantile_parameter,
                            decision_model_psa, 
                            Q=Q, 
                            N=c(5, evsi_specification$sim_data_size[i]))
    
    ### Step4(a): Simulate a future dataset from the sampling distribution p(X|phi)
    data.list.full <- list()
    
    for(j in 1:Q){
      # Generate data for moment matching
      sim_data <- data_generation_function(n=quants[j,"N"], theta=quants[j,]) 
      
      data.list.full[[j]] <- sim_data$simulated_data[decision_model_mcmc$update_prior_parameter]
    }
    
    ### Step4(b): Use Bayesian methods to update the distribution of the model 
                  # parameters
    
    # Generate Bayesian Model File
    filein <- file.path(tempdir(), fileext="datmodel.txt")
    write.model(decision_model_mcmc$model, filein)
    
    ### Step4(c): Rerun the probabilistic sensitivity analysis to update the 
                  # distribution of the net monetary benefit
    
    ### Step4(d): Calculate the variance of the net monetary benefit "sigma_q"
    
    ### Step5: Calculate the sigma^2 = Var[INB_t(theta_s)] - (1/Q)*sum(sigma_q^2)
    EVSI.var.full <- mm.post.var(model.stats=filein, data=data.list.full,
                                 N.name="n", N.size=quants[,"N"],
                                 effects=decision_model_cea$effect_funtion,
                                 costs=decision_model_cea$cost_funtion,
                                 he=bcea_result, evi=evppi.full,
                                 Q=50,
                                 data.stats=decision_model_mcmc$prior_data,
                                 update="jags",
                                 n.iter=decision_model_mcmc$iter_number,
                                 n.burnin=decision_model_mcmc$burnin_number,
                                 parameters=decision_model_mcmc$mcmc_parameter)
    
    ### Step6: Rescale the simulations of eta_t(s) such that their variance is equal to sigma^2;
    
    ### Step7: Estimates EVSI by using rescaled simulations 
    evsi.heath.full <- evsi.calc(EVSI.var.full)
    
    evsi.heath[i] <- evsi.heath.full$evsi[which(evsi.heath.full$attrib[[2]]==evsi_specification$sim_data_size[i]),
                                          which(evsi.heath.full$attrib[[1]]==decision_model_cea$wtp),
                                          which(evsi.heath.full$attrib[[3]]==0.5)]
  }
  
  return(evsi.heath)
}


