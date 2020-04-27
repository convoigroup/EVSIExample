###################################################################################
#  # A brief step-by-step guide for the Regression-Based Method #
#  1. PSA simulation to obtain theta_s and INB_t(theta_s), s = 1,...,S. 
#     (Performed by running chemotherapy_model.R).
#  2. For each s = 1,...,S:
#     (a) Simulate a dataset X_s from p(X|phi_s).
#     (b) Summarise this dataset to present how it would be analyzed in a trial,
#         denoted W(X_s).
#  3. Fit T-1 regression models with INB_t(theta_s) as the outcome and W(X_s) 
#     as the covariates.
#  4. Extract the fitted values from these regressions to estimate mu_t(X_s), 
#     which is the EVSI. 
###################################################################################

############################################################
# Function to run the Regression-Based Method
############################################################
evsi_regression_based <- function(decision_model_cea,
                                  decision_model_psa,
                                  data_generation_function,
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
  
  # Performs baseline cost-effectiveness analysis for pecified willingness to pay
  NB <- decision_model_cea$effect * decision_model_cea$wtp - decision_model_cea$cost
  
  # Find incremental net benefit for each treatment
  d.star <- which.max(apply(NB, 2, mean))
  INB <- NB - NB[, d.star]
  INB.no_reference <- INB[, -d.star, drop=F]
  
  ############################################################
  # Estimate EVSI by Strong et al method (Regression-Based Method)
  ############################################################
  
  ### Step1: PSA simulation to obtain theta_s and INB_t(theta_s)
  
  # Data for INB_t(theta_s) calculation passed to the function as 
    # decision_model_cea parameter
  # Theta_s data passed to the function as decision_model_psa parameter
  
  ### Step2(a): Simulate a dataset X_s from p(X|phi_s)
  ### Step2(b): If applicable, apply summary statistic transformation to this 
  ### dataset to present how it would be analyzed in a trial, denoted W(X_s)
  
  # Create array to store required data
  Size.Outer <- evsi_specification$outer_cond_exp_loop_size
  
  evsi.strong <- array(NA, dim=length(evsi_specification$sim_data_size))
  
  # Estimate EVSI for different sample size
  for(i in 1:length(evsi_specification$sim_data_size)) {
    print(paste0("EVSI estimation ", i, 
                 ": simulated data n=", evsi_specification$sim_data_size[i]))
    
    # Generate data
    for(j in 1:Size.Outer) {
      
      # Simulate data
      sim_data <- data_generation_function(n=evsi_specification$sim_data_size[i],
                                           theta=decision_model_psa[j,])
      
      # Convert simulated data to data frame format
      sim_data.df <- sim_data$summary_statistics[evsi_specification$gam_parameters]
      sim_data.df <- do.call(cbind, sim_data.df)
      
      # Create or append to simulated data
      if(exists("data.strong")) {
        data.strong <- rbind(data.strong, sim_data.df)
      } else {
        data.strong <- sim_data.df
      }
    }
    
    ### Step3: Fit T-1 regression models with INB_t(theta_s) as the outcome and 
               # W(X_s) as the covariates
    
    # Fit T-1 GAM regression models
    gam_formula <- paste0("INB_t ~ te(",
                          paste0(evsi_specification$gam_parameters, collapse=","),
                          ", k=", evsi_specification$gam_tensor_prod_dim, ")")
    
    # Fitted value from T-1 GAM regression models
    prepost.S <- apply(INB.no_reference[1:Size.Outer,,drop=F], c(2), 
                       function(INB_t, data.strong, gam_formula) {
                         return(gam(as.formula(gam_formula), 
                                    data=as.data.frame(data.strong))$fitted.values)
                       },
                       data.strong=data.strong, 
                       gam_formula=gam_formula)
    
    ### Step4: Extract the fitted values from these regressions to estimate mu_t(X_s), which is the EVSI
    prepost.S_max <- apply(prepost.S, c(1), max)
    
    evsi.strong[i] <- mean(pmax(0, prepost.S_max)) - max(0, mean(prepost.S_max))

    # Remove simulated data before start of next round of EVSI estimation
    rm("data.strong")
  }  
  
  return(evsi.strong)
}


