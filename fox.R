# ---------------------------------------------------------------
# JAGS code for Dynamic Biomass Production Model
# (Model with Catches and AI and 5 years projection)
# V. November 2023
# E. Rivot
# ---------------------------------------------------------------


rm(list = ls())

# ---------------------------------------------------------------

# Check the path to JAGS.exe

# Option
# If rjags can not connect directly to the JAGS.exe installed on your computer
# You will see that in the information message after loading the rjags library
# CAUTION : the path C:\\ is to be adapted to redirect where your JAGS.exe is installed
# the name JAGS-4.2.0 should also be changed if another version of JAGS is installed on your computer

# Sys.setenv(JAGS_HOME = "C:\\prog_td\\JAGS\\JAGS-4.2.0\\")

# ---------------------------------------------------------------

library(rjags)
library(coda)
library(MASS)


# Load the model
# ----------------------------------------

# Load the model, written in a .txt file
model <- "model/modelJAGS_BiomProd_PriorParamFox.txt"

# Write the model specification into a virtual text file with name "model", to be found by JAGS
# The file "model" will be used by functions "jags.model"
source(model)
model <- textConnection(model_JAGS)


# Load and format data
# -----------------------------------------

# load data

data <- read.table(file = "data/data_BiomProd.txt", header = TRUE, sep = "") 

n_obs <- length(data$Year)
n_proj <- 5
Year <- data$Year

# Vector of biomass - Needed to build the equilibrium curve
B_e <- seq(from = 0, to = 1500, by = 50)
n_equi <- length(B_e)

# Format data as a list to be read in JAGS

data <- list("I_obs" = data$I, "C_obs" = data$C, "n_obs" = n_obs, "n_proj" = n_proj,
		 "B_e" = B_e, "n_equi" = n_equi, "Year"=Year)


# MCMC options
# ----------------------------------------------------------------------

n.chains = 3

# Adapt the MCMC samplers with n.adapt iterations

n.adapt = 10000  # pre chauffe lie a JAG - avant le burning

# Iteration after adapting

n.burnin <- 1000
n.stored <- 1000
n.thin <- 10
n.iter <- n.stored*n.thin

# MANAGING NUMBER of ITERATIONS
# (works for both para = T and para = F)
# total number of REALIZED iterations = n.iter
# total number of STORED iterations = n.iter/n.thin



# Run the model to save MCMC results
# ---------------------------------------------------------------------

# Variable to store

monitor <- c( "B", "h", "D", "C", "I",
		"r", "r_p", "K", "K_p", "q", "sigma2p",
		"C_MSY", "C_MSY_p", "B_MSY", "B_MSY_p", "h_MSY", 
		"risk", "Over_C","Over_h",
		"C_e",
		"I_pred", "C_pred")

# Compile the model, create a model object with the name "model.compiled.jags" and adapt the MCMC samplers
# with n.adapt iterations

print("adapting phase")
model.compiled <- jags.model(file = model, data=data, n.chain=n.chains, n.adapt=n.adapt)

# Iteration after adapting

# Start to compute CPU time
ptm <- proc.time()

# Burnin period (not stored)

print("burn-in")
update(model.compiled, n.iter=n.burnin)

# Store mcmc samples

print("mcmc stored for results")
mcmc <- coda.samples(model=model.compiled,variable.names=monitor,n.iter=n.iter,thin=n.thin)

time.to.run.mcmc <- proc.time() - ptm
print(time.to.run.mcmc)

# ----------------------------------------------------------------------------------------------


# Run the model to compute DIC
# ---------------------------------------------------------------------

# Start from a compiled model that has already been updated

n.iter.DIC <- 1000
dic.pD <- dic.samples(model.compiled, n.iter.DIC, "pD") 
dic.pD		# Deviance Information Criterion

#Mean deviance:  202.6 
#penalty 34.93 
#Penalized deviance: 237.5 

# Alternative penalization of the Deviance
dic.popt <- dic.samples(model.compiled, n.iter, "popt")
dic.popt

#Mean deviance:  208.8 
#penalty 35.98 
#Penalized deviance: 244.8 
