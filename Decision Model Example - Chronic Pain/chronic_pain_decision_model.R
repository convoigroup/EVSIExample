# Markov Model for Chronic Pain - Adapted from Sullivan et al.
# Copyright: Anna Heath 2017
# Model has two treatment options for chronic pain and costs and utilities are 
  # determined using a Markov Model
# First-line Treatment 1: Morphine
# First-line Treatment 2: Novel Treatment
# Second-line Treatment: Oxycodone
# N: PSA sample size
# NOTE: Throughout the PSA distributions are taken as having the mean given by 
        # the parameter estimate and the 
# standard error as 10% of the parameter estimate.
# t: Number of treatment 
# l: Number of line/round


############################################################
# Load required packages
############################################################
if(!require(pacman)){               
  install.packages("pacman")
  library(pacman)
}

if(!require(BiocManager)){               
  install.packages("BiocManager")
  library(BiocManager)
}

p_install(BCEA, force=FALSE)               # Bayesian cost-effectiveness analysis
p_load(BCEA)

p_install(R2jags, force=FALSE)             # Interfaces R with 'JAGS'
p_load(R2jags)

p_install(R2OpenBUGS, force=FALSE)         # Interfaces R with 'OpenBUGS'
p_load(R2OpenBUGS)

p_install(mgcv, force=FALSE)               # Mixed GAM computation vehicle with automatic smoothness estimation
p_load(mgcv)

p_install(psych, force=FALSE)              # Procedures for psychological, psychometric, and personality research
p_load(psych)


############################################################
# Load helper functions
############################################################
##### Generate Gamma distribution parameters #####
# "par.est" is the parameter estimate, which is a prior from other literatures
# Return the value of "beta" and "scale" under "alpha" = 100

gamma.par <- function(par.est){
  alpha = 100
  beta = 100/par.est
  scale = par.est/100
  return(list(alpha = alpha,beta = beta,scale = scale))
}

##### Generate Beta distribution parameters #####
beta.par <- function(par.est){
  alpha <- (1-par.est)/0.01 - par.est
  beta <- alpha*(1/par.est - 1)
  return(list(alpha = alpha,beta = beta))
}

##### Generate Beta distribution parameters #####

betaPar <- function(m,s){
  # m - the mean
  # s - the standard deviation
  a <- m*((m*(1-m)/s^2) - 1)
  b <- (1 - m)*((m*(1-m)/s^2) - 1)
  # return beta distribution parameters
  list(a=a,b=b)
}


############################################################
# Set literature-based parameter values
############################################################
##### Parameter estimates for cycle (weekly) probabilities of tolerable AEs, 
      #treatment withdrawal and care discontinuation  #####
# First line/round of treatment: the first line treatment 
                                 # (Morphine or Novel Treatment)
# Second line/round of treatment: Oxycodone treatment if patients withdraw from 
                                  # the first treatment (due to adverse events
                                  # or lack of efficiency)
# Third line/round of treatment: further treatment if patients withdraw from the
                                 # second line treatment (due to adverse effect
                                 # or lack of efficiency)

# Cycle probability tolerable AE, treatment 1 in the first line/round (morphine)
cyc.p.ae.l1.t1 <- 0.436159243
# Cycle probability tolerable AE, treatment 2 in the first line/round (novel therapy)
cyc.p.ae.l1.t2 <- 0.324
# Cycle probability tolerable AE, treatment 1 in the second line/round (oxycodone)
cyc.p.ae.l2.t1 <- 0.4635
# Cycle probability tolerable AE, treatment 2 in the second line/round (oxycodone) 
  # (identical to treatment 1)
cyc.p.ae.l2.t2 <- cyc.p.ae.l2.t1

# The probability of experiencing adverse effect under new treatment 30% lower 
  # than that under Oxycodone
prob.novel <- (cyc.p.ae.l2.t1 - cyc.p.ae.l1.t2)/cyc.p.ae.l2.t1 
prob.novel

# Cycle probability withdraw due to AE, treatment 1 in the first line (morphine)
cyc.p.with.ae.l1.t1 <- 0.055739588 
# Cycle probability withdraw due to AE, treatment 2 in the first line (novel therapy)
cyc.p.with.ae.l1.t2 <- 0.022958454
# Cycle probability withdraw due to AE, treatment 1 in the second line (oxycodone)
cyc.p.with.ae.l2.t1 <- 0.032797792
# Cycle probability withdraw due to AE, treatment 2 in the second line (oxycodone) 
  # (identical to treatment 1)
cyc.p.with.ae.l2.t2 <- cyc.p.with.ae.l2.t1

# Cycle probability withdraw due to other reasons, treatment 1 in the first line (morphine)
cyc.p.with.l1.t1 <- 0.012741455
# Cycle probability withdraw due to other reasons, treatment 2 in the first line (novel therapy)
cyc.p.with.l1.t2 <- 0.001612408
# Cycle probability withdraw due to other reasons, treatment 1 in the second line (oxycodone)
cyc.p.with.l2.t1 <- 0.002303439
# Cycle probability withdraw due to other reasons, treatment 2 in the second line (oxycodone) 
  # (identical to treatment 1)
cyc.p.with.l2.t2 <- cyc.p.with.l2.t1

# Cycle probability of discontinuation after failed 1st-line treatment
cyc.p.discon.l1 <- 0.05
# Cycle probability of discontinuation after failed 2nd-line treatment
cyc.p.discon.l2 <- 0.1


##### Parameter estimates for model costs per cycle #####

# Treatment cost per cycle, treatment 1 in the first line/round (morphine) 
cyc.c.l1.t1 <- 2.63272916666667
# Treatment cost per cycle, treatment 2 in the first line/round (novel therapy) 
cyc.c.l1.t2 <- 55.2069

# Treatment cost per cycle, treatment 1 in the second line/round, (oxycodone)
cyc.c.l2.t1 <- 9.20115
# Treatment cost per cycle, treatment 2 in the second line/round, (oxycodone)
cyc.c.l2.t2 <- cyc.c.l2.t1

# The novel treatment is assumed to be 6 times more expensive than Oxycodone
cost.novel <- cyc.c.l1.t2/cyc.c.l2.t1
cost.novel

# Comedication cost: extra costs in order to control the adverse effect
# The comedication cost per cycle - this is for complications associated with the pain medication
# The value of comedication costs are based on a previous study, all costs have been inflated to 2013 values where appropriate
PriceIndex0910 <- 268.6
PriceIndex1213 <- 289.1
# Inflation rate
Inflation <- PriceIndex1213/PriceIndex0910 

# Co-medication cost per cycle, treatment 1 in the first line/round (morphine)
cyc.c.comed.l1.t1 <- 2.26
# Co-medication cost per cycle, treatment 2 in the first line/round (novel therapy) 
cyc.c.comed.l1.t2 <- 0.03
# Co-medication cost per cycle, treatment 1 in the second line/round, (oxycodone)
cyc.c.comed.l2.t1 <- 0.04
# Co-medication cost per cycle, treatment 2 in the second line/round, (oxycodone)
cyc.c.comed.l2.t2 <- cyc.c.comed.l2.t1

# Adverse event cost per cycle
cyc.c.ae <- 6.991009409

# The cost of withdrawing from the theraphy is identical irrespective of the reasons
# The cost of withdrawing due to adverse effect
cyc.c.withdraw.ae <- 106.911031273
# The cost of withdrawing due to other reasons (lack of efficiency)
cyc.c.withdraw. <- cyc.c.withdraw.ae

# Treatment discontinuation cost per cycle
cyc.c.discon <- 18.5


##### Parameter estimates for model utility values #####

# Utility, on treatment, no AEs (first line)
cyc.u.l1.noae <- 0.695
# Utility, on treatment, tolerable AEs (first line)
cyc.u.l1.ae <- 0.583

# Utility, withdrawn from treatment due to AEs (first line)
cyc.u.l1.withdraw.ae <- 0.503
# Utility, withdrawn from treatment due to other reasons (first line)
cyc.u.l1.withdraw.noae <- 0.405

# Utility multiplier, failed 1st-line treatment 
  # (Multiplier of the utilities under the second line treatment)
u.l2 <- 0.9
# Utility multiplier, failed 2nd-line treatment 
  # (Multiplier of the utilities under the third line treatment)
u.l3 <- 0.8


##### Time horizon #####
Time_Horizon <- 52


############################################################
# Generate PSA parameters
############################################################
##### PSA sample size #####
N <- 100000

##### Markov model transition probabilities #####

# Transition probabilities of adverse events
# Treatment 1
p.ae.l1.t1 <- rbeta(N,
                    beta.par(cyc.p.ae.l1.t1)$alpha,
                    beta.par(cyc.p.ae.l1.t1)$beta)
p.ae.l2.t1 <- rbeta(N,
                    beta.par(cyc.p.ae.l2.t1)$alpha,
                    beta.par(cyc.p.ae.l2.t1)$beta)
# Treatment 2
p.ae.l1.t2 <- p.ae.l2.t1*(1 - prob.novel) # 70% of oxycodone probability of adverse events
p.ae.l2.t2 <- p.ae.l2.t1

# Transition probilities of withdrawing due to adverse events
# Treatment 1
p.with.ae.l1.t1 <- rbeta(N,
                         beta.par(cyc.p.with.ae.l1.t1)$alpha,
                         beta.par(cyc.p.with.ae.l1.t1)$beta)
p.with.ae.l2.t1 <- rbeta(N,
                         beta.par(cyc.p.with.ae.l2.t1)$alpha,
                         beta.par(cyc.p.with.ae.l2.t1)$beta)
# Treatment 2
p.with.ae.l1.t2 <- rbeta(N,
                         beta.par(cyc.p.with.ae.l1.t2)$alpha,
                         beta.par(cyc.p.with.ae.l1.t2)$beta)
p.with.ae.l2.t2 <- p.with.ae.l2.t1

# Transition probabilities of withdrawing due to other reasons
# Treatment 1
p.with.l1.t1 <- rbeta(N,
                      beta.par(cyc.p.with.l1.t1)$alpha,
                      beta.par(cyc.p.with.l1.t1)$beta)
p.with.l2.t1 <- rbeta(N,
                      beta.par(cyc.p.with.l2.t1)$alpha,
                      beta.par(cyc.p.with.l2.t1)$beta)
# Treatment 2
p.with.l1.t2 <- rbeta(N,
                      beta.par(cyc.p.with.l1.t2)$alpha,
                      beta.par(cyc.p.with.l1.t2)$beta)
p.with.l2.t2 <- p.with.l2.t1

# Transition probability of discontinuation
p.discon.l1 <- rbeta(N,
                     beta.par(cyc.p.discon.l1)$alpha,
                     beta.par(cyc.p.discon.l1)$beta)
p.discon.l2 <- rbeta(N,
                     beta.par(cyc.p.discon.l2)$alpha,
                     beta.par(cyc.p.discon.l2)$beta)


##### Markov model transition probabilities from each state #####

# Transition probabilities for without adverse effect state for treatment 1 in the first line
                     # "No AEs" state for the 1st line\round
No.AE.l1.t1 <- cbind((1-p.with.ae.l1.t1-p.with.l1.t1)*(1-p.ae.l1.t1), 
                     # "Tolerable AEs" state for the 1st line\round
                     (1-p.with.ae.l1.t1-p.with.l1.t1)*(p.ae.l1.t1), 
                     # "Withdrawal - AEs" state for the 1st line\round
                     p.with.ae.l1.t1,  
                     # "Withdrawal - Other(including pain)" state for the 1st line\round
                     p.with.l1.t1,    
                     # "No AEs" state for the 2nd line\round
                     rep(0,N),
                     # "Tolerable AEs" state for the 2nd line\round
                     rep(0,N), 
                     # "Withdrawal - AEs" state for the 2nd line\round
                     rep(0,N), 
                     # "Withdrawal - Other(including pain)" state for the 2nd line\round
                     rep(0,N),  
                     # "Subsequent" state (absorbing state)
                     rep(0,N),  
                     # "Treatment discontinuation" state (absorbing state)
                     rep(0,N))                                        

# The sequence of these 10 states are all the same for the following transition 
  # probabilities 

# Transition probabilities for with adverse effect state for treatment 1 in the 
  # first line (Same as without adverse effect state)
AE.l1.t1 <- No.AE.l1.t1

# Transition probabilities for withdrawing due to adverse effect state for 
  # treatment 1 in the first line
With.AE.l1.t1 <- cbind(rep(0,N),
                       rep(0,N),
                       rep(0,N),
                       rep(0,N),
                       # "No AEs" state for the 2nd line\round
                       (1-p.discon.l1)*(1-p.ae.l2.t1),
                       # "Tolerable AEs" state for the 2nd line\round
                       (1-p.discon.l1)*(p.ae.l2.t1),                 
                       rep(0,N),
                       rep(0,N),
                       rep(0,N),
                       # "Treatment discontinuation" state (absorbing state)
                       p.discon.l1)                                  

# Transition probabilities for withdrawing due to other reasons (lack of efficiency) for treatment 1 in the first line (Same as withdrawing due to adverse effect state)
With.l1.t1 <- With.AE.l1.t1

# Transition probabilities for without adverse effect state for treatment 1 in the second line
No.AE.l2.t1 <- cbind(rep(0,N),
                     rep(0,N),
                     rep(0,N),
                     rep(0,N),
                     (1-p.with.ae.l2.t1-p.with.l2.t1)*(1-p.ae.l2.t1), # "No AEs" state for the 2nd line\round
                     (1-p.with.ae.l2.t1-p.with.l2.t1)*(p.ae.l2.t1),   # "Tolerable AEs" state for the 2nd line\round
                     p.with.ae.l2.t1,                                 # "Withdrawal - AEs" state for the 2nd line\round
                     p.with.l2.t1,                                    # "Withdrawal - Other(including pain)" state for the 2nd line\round
                     rep(0,N),
                     rep(0,N))
 
# Transition probabilities for with adverse effect for treatment 1 in the second line (Same as without adverse effect state)
AE.l2.t1 <- No.AE.l2.t1

# Transition probabilities for without adverse effect for treatment 2 in the first line
No.AE.l1.t2 <- cbind((1-p.with.ae.l1.t2-p.with.l1.t2)*(1-p.ae.l1.t2),   
                     (1-p.with.ae.l1.t2-p.with.l1.t2)*(p.ae.l1.t2),
                     p.with.ae.l1.t2,
                     p.with.l1.t2,
                     rep(0,N),
                     rep(0,N),
                     rep(0,N),
                     rep(0,N),
                     rep(0,N),
                     rep(0,N))

# Transition probabilities for with adverse effect for treatment 2 in the first line
AE.l1.t2 <- No.AE.l1.t2

# Transition probabilities for withdrawing due to adverse effect for treatment 2 in the first line
With.AE.l1.t2 <- cbind(rep(0,N),
                       rep(0,N),
                       rep(0,N),
                       rep(0,N),
                       (1-p.discon.l1)*(1-p.ae.l2.t2),
                       (1-p.discon.l1)*(p.ae.l2.t2),
                       rep(0,N),
                       rep(0,N),
                       rep(0,N),
                       p.discon.l1)

# Transition probabilities for withdrawing due to other reasons (lack of efficiency) for treatment 2 in the first line (Same as withdrawing due to adverse effect)
With.l1.t2 <- With.AE.l1.t2

# Transition probabilities for without adverse effect for treatment 2 in the second line
No.AE.l2.t2 <- cbind(rep(0,N),
                     rep(0,N),
                     rep(0,N),
                     rep(0,N),
                     (1-p.with.ae.l2.t2-p.with.l2.t2)*(1-p.ae.l2.t2),
                     (1-p.with.ae.l2.t2-p.with.l2.t2)*(p.ae.l2.t2),
                     p.with.ae.l2.t2,
                     p.with.l2.t2,
                     rep(0,N),
                     rep(0,N))

# Transition probabilities for with adverse effect for treatment 2 in the second line
AE.l2.t2 <- No.AE.l2.t2

# Transition probabilities for withdrawing due to adverse effect in the second line
With.AE.l2 <- cbind(rep(0,N),
                    rep(0,N),
                    rep(0,N),
                    rep(0,N),
                    rep(0,N),
                    rep(0,N),
                    rep(0,N),
                    rep(0,N),
                    rep(0,N),
                    rep(1,N))

# Transition probabilities for withdrawing due to other reasons (lack of efficiency) 
  # in the second line (Same as withdrawing due to adverse effect state)
With.l2 <- With.AE.l2

# Absorbing states: Subsequent treatment
Subs.trt <- cbind(rep(0,N),
                  rep(0,N),
                  rep(0,N),
                  rep(0,N),
                  rep(0,N),
                  rep(0,N),
                  rep(0,N),
                  rep(0,N),
                  rep(1,N),
                  rep(0,N))

# Absorbing states: Treatment discontinuation
trt.discon <- With.AE.l2


##### Markov model transition matrices #####

# Setup a 3-dimensional array containing transition matrix for each PSA simulation
PSA.Trans.Mat.t1 <- array(NA,dim=c(10,10,N))

# Building transition matrix by combining transition probabilities from all states
# Treatment 1
  # Fill the first dimension of the transition matrix with a 10 x N matrix for 10 states
  # Use transpose since we need a 10 x N matrix but transition probabilities are all N x 10

# "No AEs" state for the 1st line\round
PSA.Trans.Mat.t1[1,,] <- (t(No.AE.l1.t1)) 
# "Tolerable AEs" state for the 1st line\round
PSA.Trans.Mat.t1[2,,] <- (t(AE.l1.t1)) 
# "Withdrawal - AEs" state for the 1st line\round
PSA.Trans.Mat.t1[3,,] <- (t(With.AE.l1.t1)) 
# "Withdrawal - Other(including pain)" state for the 1st line\round
PSA.Trans.Mat.t1[4,,] <- (t(With.l1.t1)) 
# "No AEs" state for the 2nd line\round
PSA.Trans.Mat.t1[5,,] <- (t(No.AE.l2.t1))  
# "Tolerable AEs" state for the 2nd line\round
PSA.Trans.Mat.t1[6,,] <- (t(AE.l2.t1)) 
# "Withdrawal - AEs" state for the 2nd line\round
PSA.Trans.Mat.t1[7,,] <- (t(With.AE.l2)) 
# "Withdrawal - Other(including pain)" state for the 2nd line\round
PSA.Trans.Mat.t1[8,,] <- (t(With.l2))   
# "Subsequent" state (absorbing state)
PSA.Trans.Mat.t1[9,,] <- (t(Subs.trt)) 
# "Treatment discontinuation" state (absorbing state)
PSA.Trans.Mat.t1[10,,] <- (t(trt.discon))         

# Treatment 2
PSA.Trans.Mat.t2 <- PSA.Trans.Mat.t1
PSA.Trans.Mat.t2[1,,] <- (t(No.AE.l1.t2))
PSA.Trans.Mat.t2[2,,] <- (t(AE.l1.t2))
PSA.Trans.Mat.t2[3,,] <- (t(With.AE.l1.t2))
PSA.Trans.Mat.t2[4,,] <- (t(With.l1.t2))
PSA.Trans.Mat.t2[5,,] <- (t(No.AE.l2.t2))
PSA.Trans.Mat.t2[6,,] <- (t(AE.l2.t2))


############################################################
# Run model simulation
############################################################
##### State-Transition Markov cohort model #####
Markov_Prob <- function(TransMat) {
  # Simulation without half cycle correction
  # Set up a 52 x 10 trace matrix containing probabilities of each state at 
    # each time point
  trace_matrix <- matrix(NA, nrow = Time_Horizon, ncol = ncol(TransMat))
  
  # Initial state 
    # Everyone starts from the first state
  InitVector <- c(1,0,0,0,0,0,0,0,0,0)
  
  # Set and simulate state probabilities for the first two time points
  # The state probabilities of current time point equals the state probabilities
    # of previous time point multiply by "TransMat"
  trace_matrix[1,] <- InitVector
  
  # Simulate state probabilities for the remaining time points
  for (i in 2:nrow(trace_matrix)) {
    trace_matrix[i,] <- trace_matrix[i-1,] %*% TransMat
  }
  
  # Simulation with half cycle correction (hcc)
    # Due to the short time horizons considered, model health outcomes and costs
      # were not discounted, so a half-cycle correction was applied 
  
  # Setup matrix containing state probabilities at each time point
  hcc_trace_matrix <- matrix(NA, nrow = Time_Horizon, ncol = ncol(TransMat))
  
  # Apply hcc correction to state probabilities for the first time point
  hcc_trace_matrix[1,] <- 0.5*InitVector + 0.5*trace_matrix[2,]
  
  # Apply hcc correction to state probabilities for the remaining time point
  for (i in 2:Time_Horizon-1) {
    hcc_trace_matrix[i,] <- 0.5*trace_matrix[i,] + 0.5*trace_matrix[i+1,]
  }
  
  # No hcc correction to state probabilities for the last time point
  hcc_trace_matrix[Time_Horizon,] <- trace_matrix[Time_Horizon,]
  
  # Output hcc-corrected state probabilities
  return(hcc_trace_matrix)
}


##### Markov cohort model simulation results #####
# "Prob.Array.1" stores the Markov simulation results for each PSA simulation,
  # it's an array of 3-dimensions: 
    # time points(52) x states(10) x PSA simulations(N)

# Probabilities for each state at each time point for treatment 1
Prob.Array.1 <- array(NA,c(Time_Horizon,10,N))
for(i in 1:N){
  Prob.Array.1[,,i] <- Markov_Prob(PSA.Trans.Mat.t1[,,i])
}

# Probabilities for each state at each time point for treatment 2
Prob.Array.2 <- array(NA,c(Time_Horizon,10,N))
for(i in 1:N){
  Prob.Array.2[,,i] <- Markov_Prob(PSA.Trans.Mat.t2[,,i])
}


############################################################
# Estimate cost
############################################################
##### Generate PSA parameters for cost #####
# Treatment costs are considered known and the values were taken from the literature

# Cost of treatment 1 (Morphine) of the first round
c.l1.t1 <- cyc.c.l1.t1 
# Cost of treatment 2 (Novel Treatment) of the second round
c.l1.t2 <- cyc.c.l1.t2

# Comedication cost of treatment 1 (Morphine)
c.comed.l1.t1 <- rgamma(N,
                        shape = gamma.par(cyc.c.comed.l1.t1*Inflation)$alpha,
                        rate = gamma.par(cyc.c.comed.l1.t1*Inflation)$beta)

# Comedication cost of treatment 2 (Novel Treatment)
  # -- is reduced by 30% compared to oxycodone
c.comed.l1.t2 <- rgamma(N,
                        shape = gamma.par(cyc.c.comed.l2.t1*Inflation*(1-prob.novel))$alpha,
                        rate = gamma.par(cyc.c.comed.l2.t1*Inflation*(1-prob.novel))$beta)

# Costs of adverse events
c.ae <- rgamma(N,
               shape = gamma.par(cyc.c.ae)$alpha,
               rate = gamma.par(cyc.c.ae)$beta)

# The cost of withdrawing due to adverse effect
c.withdraw.ae <- rgamma(N,
                        shape = gamma.par(cyc.c.withdraw.ae)$alpha,
                        rate = gamma.par(cyc.c.withdraw.ae)$beta)

# The cost of withdrawing due to lack of efficiency
# The cost of withdrawing from the theraphy is the same irrespective of the reasons
c.withdraw <- c.withdraw.ae

# Cost of discontinuing treatment is based on visiting the GP (general
  # practitioner surgery visit)
c.discon <- rgamma(N,
                   shape = gamma.par(cyc.c.discon)$alpha,
                   rate = gamma.par(cyc.c.discon)$beta)

# Cost of the second round of the treatment 
  # (adjusted by comedication cost and the inflation rate)
c.l2 <- cyc.c.l2.t1 + cyc.c.comed.l2.t1 * Inflation

# Cost of the third round of the treatment 
  # (adjusted by comedication cost and the inflation rate)
c.l3 <- cyc.c.l1.t1 + cyc.c.comed.l1.t1 * Inflation


##### Cost matrices containing all of the costs in all three lines/rounds #####
# Cost matrices (N x 10) for each state in treatment 1
c.Mat.t1 <- cbind(c.l1.t1 + c.comed.l1.t1,
                  c.l1.t1 + c.comed.l1.t1 + c.ae,
                  c.withdraw.ae,
                  c.withdraw,
                  c.l2,
                  c.l2 + c.ae,
                  c.withdraw.ae,
                  c.withdraw,
                  c.l3,
                  c.discon)

# Cost matrices (N x 10) for each state in treatment 2
c.Mat.t2 <- cbind(cyc.c.l1.t2 + c.comed.l1.t2,
                  cyc.c.l1.t2 + c.comed.l1.t2 + c.ae, 
                  c.withdraw.ae,
                  c.withdraw,
                  c.l2,
                  c.l2 + c.ae,
                  c.withdraw.ae,
                  c.withdraw,
                  c.l3,
                  c.discon)


##### Cost estimation #####

# Total cost of 52 time points for each PSA simulation for treatment 1
costs.t1 <- array(NA,N)
for(i in 1:N){
  costs.t1[i] <- sum(Prob.Array.1[,,i] %*% c.Mat.t1[i,])
}

# Total cost of of 52 time points for each PSA simulation for treatment 2
costs.t2 <- array(NA,N)
for(i in 1:N){
  costs.t2[i] <- sum(Prob.Array.2[,,i] %*% c.Mat.t2[i,])
}


############################################################
# Estimate utility
############################################################

##### Generate PSA parameters for utilities #####

# The utilities of the first line treatment without adverse events 
u.l1.noae <- rbeta(N,
                   beta.par(cyc.u.l1.noae)$alpha,
                   beta.par(cyc.u.l1.noae)$beta)

# The utilities of the first line treatment with adverse events 
u.l1.ae <- rbeta(N,
                 beta.par(cyc.u.l1.ae)$alpha,
                 beta.par(cyc.u.l1.ae)$beta)

# The utilities of withdrawing from the first line treatment due to adverse events
u.l1.withdraw.ae <- rbeta(N,
                          beta.par(cyc.u.l1.withdraw.ae)$alpha,
                          beta.par(cyc.u.l1.withdraw.ae)$beta)

# The utilities of withdrawing from the first line treatment due to other reasons
u.l1.withdraw.noae <- rbeta(N,
                            beta.par(cyc.u.l1.withdraw.noae)$alpha,
                            beta.par(cyc.u.l1.withdraw.noae)$beta)

# Multiplier of the utilities under the second line treatment 
u.l2 <- u.l2

# The utility for the third line treatment
u.l3 <- (u.l1.noae + u.l1.ae)/2

# The utility for discontinuing the treatment 
# Treatment discontinuation is attributed the lowest utility estimate (utility 
  # of withdrawn from treatment due to other reasons), adjusted by an assumed 
  # multiplier of 0.8 to reflect the negative effect of successive failed treatments 
  # upon HRQL (health-related quality of life)
u.discon <- u.l1.withdraw.noae*0.8


##### Utility matrices containing all of the utilities in all three lines per week #####
u.Mat <- cbind(u.l1.noae,
               u.l1.ae,
               u.l1.withdraw.ae,
               u.l1.withdraw.noae,
               u.l1.noae*u.l2,
               u.l1.ae*u.l2,
               u.l1.withdraw.ae*u.l2,
               u.l1.withdraw.noae*u.l2,
               u.l3,
               u.discon)*7/365.25


##### Utility estimation #####

# Total utility of 52 time points for each PSA simulation for treatment 1
utility.t1 <- array(NA,N)
for(i in 1:N){
  utility.t1[i] <- sum(Prob.Array.1[,,i]%*%u.Mat[i,])
}

# Total utility of 52 time points for each PSA simulation for treatment 2
utility.t2 <- array(NA,N)
for(i in 1:N){
  utility.t2[i] <- sum(Prob.Array.2[,,i]%*%u.Mat[i,])
}


############################################################
# Perform baseline cost and effectiveness analysis
############################################################

discount.15 <- sum(1/(1+0.035)^(0:15))

m <- bcea(discount.15*cbind(utility.t1,utility.t2),
          discount.15*cbind(costs.t1,costs.t2),
          ref=2,
          wtp=c(0,20000))

# Prediction variance for a specific willingness to pay
var.pr <- var(m$ib[which(m$k==20000),])


############################################################
# Package cost, utility and model parameters for output
############################################################
pars <- cbind(c.ae,
              c.discon,
              c.comed.l1.t1,
              c.comed.l1.t2,
              c.withdraw,
              c.withdraw.ae,
              p.ae.l1.t1,
              p.ae.l1.t2,
              p.discon.l1,
              p.discon.l2,
              p.with.ae.l1.t1,
              p.with.ae.l1.t2,
              p.with.l1.t1,
              p.with.l1.t2,
              p.ae.l2.t1,
              p.ae.l2.t2,
              p.with.ae.l2.t1,
              p.with.ae.l2.t2,
              p.with.l2.t1,
              p.with.l2.t2,
              u.discon,
              u.l1.ae,
              u.l1.noae,
              u.l1.withdraw.ae,
              u.l1.withdraw.noae,
              u.l3)

