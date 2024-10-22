library(here)
library(nimble)
library(rcarbon)
library(nimbleCarbon)
source(here('src','utility.R'))
data(intcal20)
seed  <- 123

# Define Simulation model ----
diffusionModel  <- nimbleCode({
	for (i in 1:N)
	{
		y[i] ~ dbern(p[i]) #probability of event y
		p[i]  <- k[siteID[i]]/(1+exp(r*(theta[i]-m))); #sigmoidal model


		#Calibration of observed cra
		mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigma[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=mu[i],sd=sigma[i]);

		# Prior for theta (flat)
		theta[i] ~ dunif(1000,10000)
	}

	r ~ dexp(10) # prior diffusion rate
	m ~ T(dnorm(mean=midPrior,sd=midSD),1000,50000) #prior mid-point
	for (j in 1:NSites)
	{
		k[j] ~ dbeta(beta0,beta1)
	}
	beta0  <- mu_k * phi + 1
	beta1  <- (1-mu_k) * phi + 1
})


# Simulation 1a (based on Case Study 1a - Japan) ----
# Model Parameters
r  <- 0.01
m  <- 2900
mu_k  <- 0.65
phi  <- 50

true.param.1a  <- list(r=r,m=m,mu_k=mu_k,phi=phi)

# Prior Settings for m
midPrior  <- 3200
midSD  <- 500

# Sample sizes
nsites  <- 215
nsamples  <- 551

# Time window
begin  <- 4000
end  <- 1700

# Define Constants
constants  <- list()
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma
constants$midPrior  <- midPrior
constants$midSD  <- midSD
constants$N  <- nsamples
constants$NSites  <- nsites
constants$siteID  <- c(1:constants$NSites,sample(1:constants$NSites,size=constants$N-constants$NSites,replace=TRUE)) #Assign SiteIDs
constants$cra_error  <- rep(20,constants$N) #Constant 14C error of 20yrs
constants$mu_k  <- mu_k
constants$phi  <- phi

# Simulate Response Variable
diffusionModel.sim  <- nimbleModel(diffusionModel,constants=constants)
set.seed(seed)
diffusionModel.sim$theta  <- round(runif(constants$N,min=end,max=begin)) #Random true dates uniformly within time-window
diffusionModel.sim$r  <- r
diffusionModel.sim$m  <- m
diffusionModel.sim$simulate('mu')
diffusionModel.sim$simulate('sigmaCurve')
diffusionModel.sim$simulate('sigma')
diffusionModel.sim$simulate('cra')
diffusionModel.sim$calculate('beta0')
diffusionModel.sim$calculate('beta1')
diffusionModel.sim$simulate('k')
diffusionModel.sim$simulate('p')
diffusionModel.sim$simulate('y')

# Store Output
d  <- list()
d$y  <- diffusionModel.sim$y
d$cra  <- round(diffusionModel.sim$cra)
save(constants,d,true.param.1a,file=here('sim','simdata','simdata1a.RData'))


# Simulation 1b (based on Case Study 1b - Britain) ----
# Model Parameters
r  <- 0.008
m  <- 4500
mu_k  <- 0.8
phi  <- 40
true.param.1b  <- list(r=r,m=m,mu_k=mu_k,phi=phi)
# Prior Settings for m
midPrior  <- 3500
midSD  <- 500

# Sample sizes
nsites  <- 314
nsamples  <- 928

# Time window
begin  <- 7000
end  <- 3000

# Define Constants
constants  <- list()
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma
constants$midPrior  <- midPrior
constants$midSD  <- midSD
constants$N  <- nsamples
constants$NSites  <- nsites
constants$siteID  <- c(1:constants$NSites,sample(1:constants$NSites,size=constants$N-constants$NSites,replace=TRUE)) #Assign SiteIDs
constants$cra_error  <- rep(20,constants$N) #Constant 14C error of 20yrs
constants$mu_k  <- mu_k
constants$phi  <- phi

# Simulate Response Variable
diffusionModel.sim  <- nimbleModel(diffusionModel,constants=constants)
set.seed(seed)
diffusionModel.sim$theta  <- round(runif(constants$N,min=end,max=begin)) #Random true dates uniformly within time-window
diffusionModel.sim$r  <- r
diffusionModel.sim$m  <- m
diffusionModel.sim$simulate('mu')
diffusionModel.sim$simulate('sigmaCurve')
diffusionModel.sim$simulate('sigma')
diffusionModel.sim$simulate('cra')
diffusionModel.sim$calculate('beta0')
diffusionModel.sim$calculate('beta1')
diffusionModel.sim$simulate('k')
diffusionModel.sim$simulate('p')
diffusionModel.sim$simulate('y')

# Store Output
d  <- list()
d$y  <- diffusionModel.sim$y
d$cra  <- round(diffusionModel.sim$cra)
save(constants,d,true.param.1b,file=here('sim','simdata','simdata1b.RData'))
