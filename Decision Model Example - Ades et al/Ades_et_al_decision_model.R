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

p_install(boot, force=FALSE)               # Bootstrap functions
p_load(boot)


############################################################
# Load helper functions
############################################################


############################################################
# Set literature-based parameter values
############################################################
L <- 30
Qse <- 1
Ce <- 200000
Ct <- 15000
Cse <- 100000
lambda <- 75000


############################################################
# Generate PSA parameters
############################################################

##### PSA sample size #####
N <- 1000000

lQe <- rnorm(N,0.6,sqrt(1/6))
Qe <- inv.logit(lQe)
Pc <- rbeta(N,15,85)
Pse <- rbeta(N,3,9)
OR <- exp(rnorm(N,-1.5,sqrt(1/3)))
Pt <- inv.logit(log(Pc/(1-Pc))+log(OR))


############################################################
# Run model simulation
############################################################


############################################################
# Estimate cost
############################################################

cost.standard <- 
  Pc*(Ce)+
  (1-Pc)*0

cost.treatment <- 
  Pse*Pt*((Ct+Cse+Ce))+
  Pse*(1-Pt)*((Ct+Cse))+
  (1-Pse)*Pt*((Ct+Ce))+
  (1-Pse)*(1-Pt)*(Ct)


############################################################
# Estimate utility
############################################################

utility.standard <- 
  Pc*(L*(1+Qe)/2)+
  (1-Pc)*L

utility.treatment <- 
  Pse*Pt*((L*(1+Qe)/2-Qse))+
  Pse*(1-Pt)*((L-Qse))+
  (1-Pse)*Pt*(L*(1+Qe)/2)+
  (1-Pse)*(1-Pt)*(L)


############################################################
# Perform baseline cost and effectiveness analysis
############################################################


############################################################
# Package cost, utility and model parameters for output
############################################################
cost <- data.frame("treatment"=cost.treatment,
                   "standard"=cost.standard)

utility <- data.frame("treatment"=utility.treatment,
                      "standard"=utility.standard)

theta <- data.frame(L=L,
                    lQe=lQe,
                    Qe=Qe,
                    Qse=Qse,
                    Ce=Ce,
                    Ct=Ct,
                    Cse=Cse,
                    Pc=Pc,
                    Pse=Pse,
                    OR=OR,
                    Pt=Pt,
                    lambda=lambda)


