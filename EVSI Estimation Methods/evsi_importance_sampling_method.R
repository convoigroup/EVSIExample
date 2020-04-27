###################################################################################
#  # A brief step-by-step guide for the Importance Sampling Method #
#  1. PSA simulation to obtain theta_s and INB_t(theta_s), s = 1,...,S.
#     (Performed by running chemotherapy_model.R).
#  2. Estimate EVPPI to obtain simulations of the inner conditional expectation
#     # of INB on the future data (eta_t(s)).
#  3. For each s = 1,...,S:
#     (a) Simulate a dataset X_s from p(X|phi_s).
#     (b) Compute the likelihood of X_s conditional on all PSA simulations for 
#         theta, L_r, r = 1,...,S.
#     (c) Compute l_r = L_r/[Sum(L_r)], so l_r sums to 1.
#     (d) Calculate the weighted sum of eta_t(r), Sum[l_r * eta_t(r)].
#  4. Each weighted sum estimates mu_t(X_s), which is the EVSI. 
###################################################################################

############################################################
# Function to run the Importance Sampling Method
############################################################
evsi_importance_sample <- function(decision_model_cea,
                                   decision_model_psa,
                                   data_generation_function,
                                   likelihood_compute_function,
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
  
  # Performs baseline cost-effectiveness analysis for specified willingness to pay
  NB <- decision_model_cea$effect*decision_model_cea$wtp - decision_model_cea$cost
  
  # Find incremental net benefit for each treatment
  d.star <- which.max(apply(NB, 2, mean))
  INB <- NB - NB[, d.star]
  INB.no_reference <- INB[, -d.star, drop=F]
  
  ############################################################
  # Estimate EVSI by the Importance Sampling Method
  ############################################################
  
  ### Step1: PSA simulation to obtain theta_s and INB_t(theta_s)
  
  ### Step2: Estimate EVPPI to obtain simulations of the inner conditional 
             # expectation of INB on the future data (eta_t(s))

  # Estimate EVPPI to obtain the inner conditional expectation of INB on the future data
  gam_formula <- paste0("INB_t ~ te(",
                        paste0(evsi_specification$evppi_parameters, collapse=","),
                        ", k=", evsi_specification$gam_tensor_prod_dim, ")")
  
  # Conditional distribution on the parameter
  con.dist <- apply(INB.no_reference, c(2), 
                     function(INB_t, decision_model_psa, gam_formula) {
                       return(gam(as.formula(gam_formula), 
                                  data = decision_model_psa)$fitted.values)
                     },
                    decision_model_psa = decision_model_psa, 
                    gam_formula = gam_formula)

  ### Step3
  
  # Create array to store required data
  Size.Outer <- evsi_specification$outer_cond_exp_loop_size

  prepost.M <- array(dim=Size.Outer)
  
  evsi.menzies <- array(NA, dim=length(evsi_specification$sim_data_size))
  
  # Calculate EVSI for different sample size
  for(i in 1:length(evsi_specification$sim_data_size)) {
    print(paste0("EVSI estimation ", i, 
                 ": simulated data n=", evsi_specification$sim_data_size[i]))
    
    for(j in 1:Size.Outer){

      ### Step3(a): Simulate a dataset X_s from p(X|phi_s)
      # Simulate dataset
      sim_data <- data_generation_function(n=evsi_specification$sim_data_size[i],
                                           theta=decision_model_psa[j,]) 
      
      # Convert simulated data to data frame format
      sim_data.df <- sim_data$simulated_data[evsi_specification$likelihood_data_variables]
      sim_data.df <- do.call(cbind, sim_data.df)
      
      ### Step3(b): Compute the likelihood of X_s conditional on all PSA 
                    # simulations for theta, L_r, r = 1,...,S
      likelihood <- likelihood_compute_function(n=evsi_specification$sim_data_size[i],
                                                D=sim_data.df,
                                                theta=decision_model_psa[1:Size.Outer,]) 

      ### Step3(c): Compute l_r = L_r/[Sum(L_r)], so l_r sums to 1
      weight <- likelihood / sum(likelihood)               
      
      ### Step3(d): Calculate the weighted sum of eta_t(r), Sum[l_r * eta_t(r)]
      prepost.M[j] <- weight %*% con.dist[1:Size.Outer,]
    }
    
    ### Step4: Each weighted sum estimates mu_t(X_s), which is the EVSI.
    prepost.M_max <- apply(prepost.M, c(1), max)
    
    evsi.menzies[i] <- mean(pmax(0, prepost.M_max), na.rm=TRUE) - 
      max(0, mean(prepost.M_max), na.rm=TRUE)
  }
  
  return(evsi.menzies)
}


