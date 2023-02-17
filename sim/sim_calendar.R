library(here)
library(nimble)
library(nimbleCarbon)
source(here('src','utility.R'))


# Define model ----
adoptionModel  <- nimbleCode({
	for (i in 1:N)
	{
		y[i] ~ dbern(p[i]) #probability of event y
		p[i]  <- k[siteID[i]]/(1+exp(r*(BP[i]-m))); #sigmoidal model
	}

	r ~ dexp(10) # prior adoption rate
	m ~ T(dnorm(mean=midPrior,sd=midSD),1000,50000) #prior mid-point
	for (j in 1:NSites)
	{
		logk[j] ~ dnorm(mean=mu_k,sd=sigma_k) #prior site
		k[j]  <- 1/(1+exp(-logk[j])) 
	}
	mu_k ~ dnorm(0,1) #hyperprior for site prior
	sigma_k ~ dinvgamma(5,5) #hyperprior for site prior
})

# Prior Predictive checks ----
# k
hist(logistic(rnorm(1000,mean=rnorm(1000,0,1),sd=rinvgamma(1000,5,5))))


# Model Parameters ----
r  <- 0.004
m  <- 2700
mu_k  <- 1
sigma_k  <- 0.5


# Simulate ----
# Define Constants
constants  <- list()
constants$midPrior  <- 3000 
constants$midSD  <- 500
constants$N  <- 400
constants$NSites  <- 100
constants$siteID  <- c(1:constants$NSites,sample(1:constants$NSites,size=constants$N-constants$NSites,replace=TRUE))
constants$BP  <- round(runif(constants$N,max=5000,1500))

# Simulate Response Variable
adoptionModel.sim  <- nimbleModel(adoptionModel,constants=constants)
adoptionModel.sim$r  <- r
adoptionModel.sim$m  <- m
adoptionModel.sim$mu_k  <- m_k
adoptionModel.sim$sigma_k  <- sigma_k
adoptionModel.sim$simulate('logk')
adoptionModel.sim$simulate('k')
adoptionModel.sim$simulate('p')
# plot(sort(constants$BP),adoptionModel.sim$p[order(constants$BP)],xlim=c(5000,1500))
adoptionModel.sim$simulate('y')

# Store simulated variables 
d  <- list()
d$y  <- adoptionModel.sim$y

# Fit ----
inits  <- list()
inits$r  <- 0.0001
inits$m  <- 3000
inits$mu_k  <- 1
inits$sigma_k  <- 1
inits$logk  <- rnorm(constants$NSites,mean=inits$mu_k,sd=inits$sigma_k)
fit.model  <- nimbleModel(adoptionModel,constants=constants,data=d,inits=inits)
cfit.model <- compileNimble(fit.model)
conf <- configureMCMC(fit.model)
MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC)
results <- runMCMC(cMCMC, niter = 200000,nchains=3, thin=10,nburnin = 100000,samplesAsCodaMCMC = T)

# Diagnostic and Posterior Plots ----
coda::gelman.diag(results)
res  <- do.call(rbind.data.frame,results)
par(mfrow=c(2,2))
postHPDplot(res[,'r'])
abline(v=r)
postHPDplot(res[,'m'])
abline(v=m)
postHPDplot(res[,'mu_k'])
abline(v=mu_k)
postHPDplot(res[,'sigma_k'])
abline(v=sigma_k)
