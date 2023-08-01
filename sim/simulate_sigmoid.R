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
# 		logk[j] ~ dnorm(mean=mu_k,sd=sigma_k) #prior site
# 		k[j]  <- 1/(1+exp(-logk[j])) 
		k[j] ~ dbeta(beta0,beta1)
	}
# 	mu_k ~ dnorm(0,1) #hyperprior for site prior
# 	sigma_k ~ dexp(10)
# 	mu_k ~ dbeta(2,2) #Unimodal shape with low prob on 0 and 1
# 	sigma_k ~ dnorm(0,0.0001)
# 	tau  <- 1/sqrt(sigma_k) 
	beta0  <- mu_k * tau + 1
	beta1  <- (1-mu_k) * tau + 1
# 	sigma_k ~ dinvgamma(5,5) #hyperprior for site prior
})


# Simulation 1 (based on Case Study I - Japan) ----
# Model Parameters
r  <- 0.004
m  <- 2700
mu_k  <- 0.8
tau  <- 20

true.param.1  <- list(r=r,m=m,mu_k=mu_k,tau=tau)

# Prior Settings for m
midPrior  <- 3000
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
constants$tau  <- tau

# Simulate Response Variable
diffusionModel.sim  <- nimbleModel(diffusionModel,constants=constants)
set.seed(seed)
diffusionModel.sim$theta  <- round(runif(constants$N,min=end,max=begin)) #Random true dates uniformly within time-window
diffusionModel.sim$r  <- r
diffusionModel.sim$m  <- m
# diffusionModel.sim$mu_k  <- mu_k
# diffusionModel.sim$sigma_k  <- sigma_k
diffusionModel.sim$simulate('mu')
diffusionModel.sim$simulate('sigmaCurve')
diffusionModel.sim$simulate('sigma')
diffusionModel.sim$simulate('cra')
# diffusionModel.sim$simulate('logk')
diffusionModel.sim$calculate('beta0')
diffusionModel.sim$calculate('beta1')
diffusionModel.sim$simulate('k')
diffusionModel.sim$simulate('p')
plot(sort(diffusionModel.sim$theta),diffusionModel.sim$p[order(diffusionModel.sim$theta)],xlim=c(5000,1500))
diffusionModel.sim$simulate('y')

# Store Output
d  <- list()
d$y  <- diffusionModel.sim$y
d$cra  <- round(diffusionModel.sim$cra)
save(constants,d,true.param.1,file=here('sim','simdata','simdata1.RData'))


# Simulation 2 (based on Case Study I - Britain) ----
# Model Parameters
r  <- 0.008
m  <- 4500
mu_k  <- 0.9
tau  <- 40
true.param.2  <- list(r=r,m=m,mu_k=mu_k,tau=tau)
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
constants$tau  <- tau

# Simulate Response Variable
diffusionModel.sim  <- nimbleModel(diffusionModel,constants=constants)
set.seed(seed)
diffusionModel.sim$theta  <- round(runif(constants$N,min=end,max=begin)) #Random true dates uniformly within time-window
diffusionModel.sim$r  <- r
diffusionModel.sim$m  <- m
# diffusionModel.sim$mu_k  <- mu_k
# diffusionModel.sim$sigma_k  <- sigma_k
diffusionModel.sim$simulate('mu')
diffusionModel.sim$simulate('sigmaCurve')
diffusionModel.sim$simulate('sigma')
diffusionModel.sim$simulate('cra')
diffusionModel.sim$calculate('beta0')
diffusionModel.sim$calculate('beta1')
# diffusionModel.sim$simulate('logk')
diffusionModel.sim$simulate('k')
diffusionModel.sim$simulate('p')
plot(sort(diffusionModel.sim$theta),diffusionModel.sim$p[order(diffusionModel.sim$theta)],xlim=c(7000,3000))
diffusionModel.sim$simulate('y')

# Store Output
d  <- list()
d$y  <- diffusionModel.sim$y
d$cra  <- round(diffusionModel.sim$cra)
save(constants,d,true.param.2,file=here('sim','simdata','simdata2.RData'))
