# Load libraries and Data ----
library(here)
library(parallel)
library(nimbleCarbon)
library(rcarbon)
data(intcal20)
load(here('data','jpdata.RData'))

# Data Preparation ----
caldates  <- calibrate(jpdata$C14Age,jpdata$C14Error)
caldates  <- caldates[ii]
jpdata$site.id  <- as.integer(factor(jpdata$SiteNameJp))

# Define Constants ----
constants  <- list()
constants$N  <- nrow(jpdata)
constants$NSites  <- length(unique(jpdata$site.id))
constants$siteID  <- jpdata$site.id
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma
constants$cra_error  <- jpdata$C14Error

# Define Data ----
d  <- list()
d$y  <- jpdata$rice_nuts 
d$cra  <- jpdata$C14Age

# Define Inits for theta ----
theta  <- medCal(caldates)

# Core Runscript  ----
runFun  <- function(seed, d, constants, theta, init, nburnin, niter, thin)
{
	#Load R libraries
	library(nimbleCarbon)
	library(truncnorm)

	#Adoption model
	adoptionModel  <- nimbleCode({
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
			theta[i] ~ dunif(1000,50000)
		}

		r ~ dexp(10) # prior adoption rate
		m ~ T(dnorm(mean=2500,sd=500),1000,50000) #prior mid-point
		for (j in 1:NSites)
		{
			logk[j] ~ dnorm(mean=mu_k,sd=sigma_k) #prior site
			k[j]  <- 1/(1+exp(-logk[j])) 
		}
		mu_k ~ dnorm(0,1) #hyperprior for site prior
		sigma_k ~ dinvgamma(5,5) #hyperprior for site prior
	})

	#Define inits
	inits  <- list()
	inits$r  <- 0.0001
	inits$m  <- 3000
	inits$mu_k  <- 1
	inits$sigma_k  <- 1
	inits$logk  <- rnorm(constants$NSites,mean=inits$mu_k,sd=inits$sigma_k)
	inits$theta  <- theta

	#Setup MCMC
	fit.model  <- nimbleModel(adoptionModel,constants=constants,data=d,inits=inits)
	cfit.model <- compileNimble(fit.model)
	conf <- configureMCMC(fit.model)
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)

	results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed)
}

# MCMC Setup ----
niter  <- 1000000
nchains  <- 4
nburnin  <- niter/2
thin  <- ceiling((niter-nburnin)/10000)
cl  <- makeCluster(nchains)
seeds  <- c(12,34,56,78)[1:nchains]

# Run MCMC ----
out  <- parLapply(cl=cl,X=seeds,fun=runFun,d=d,constants=constants,theta=theta,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(cl)

# Diagnostic and Posterior Processing ----
post.sample  <- coda::mcmc.list(out)
rhats.jp.abot  <- coda::gelman.diag(post.sample)
# which(rhat.jp.abot[[1]][,1]>1.01)

# Store output ----
post.sample.combined  <- do.call(rbind.data.frame,post.sample)
post.sample.theta  <- post.sample.combined[,grep('theta',colnames(post.sample.combined))]
post.sample.core.jp.abot  <- post.sample.combined[,!grepl('theta',colnames(post.sample.combined))]


# Posterior Predictive Checks ----
nsim  <- 1000 #Number of posterior simulations
ppmat  <- matrix(NA,nrow=constants$N,ncol=nsim) #Matrix storing predictions

# Simulation Model
adoptionSimModel  <- nimbleCode({
	for (i in 1:N)
	{
		y[i] ~ dbern(p[i])
		p[i]  <- k[siteID[i]]/(1+exp(r*(theta[i]-m))); #sigmoidal model
		mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigma[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=mu[i],sd=sigma[i]);
		theta[i] ~ dunif(1000,10000)
	}
	r ~ dexp(100) # prior adoption rate
	m ~ T(dnorm(mean=5500,sd=1000),1000,50000) #prior mid-point

	for (j in 1:NSites)
	{
		logk[j] ~ dnorm(mean=mu_k,sd=sigma_k) #prior site
		k[j]  <- 1/(1+exp(-logk[j])) 
	}
	mu_k ~ dnorm(0,1) #hyperprior for site prior for region 1
	sigma_k ~ dinvgamma(5,5) #hyperprior for site prior
})

s.index  <- sample((niter-nburnin)*nchains/thin,size=nsim)

sim.model  <- nimbleModel(adoptionSimModel,constants=constants,data=d)

pb <- txtProgressBar(min = 0, max = nsim, style = 3, width = 50, char = "=")

for (i in 1:nsim)
{
    setTxtProgressBar(pb, i)
    ii  <- s.index[i]
    sim.model$sigma_k  <- as.numeric(post.sample.combined[ii,'sigma_k'])
    sim.model$mu_k  <- as.numeric(post.sample.combined[ii,'mu_k'])
    sim.model$r  <- as.numeric(post.sample.combined[ii,'r'])
    sim.model$m  <- as.numeric(post.sample.combined[ii,'m'])
    sim.model$theta  <- as.numeric(post.sample.combined[ii,grep('theta\\[',colnames(post.sample.combined))])
    sim.model$simulate('logk')
#     sim.model$logk  <- as.numeric(post.sample.combined[ii,grep('logk\\[',colnames(post.sample.combined))])
    sim.model$calculate('k')
    sim.model$simulate('p')
    ppmat[,i]  <- rbinom(constants$N,prob=unlist(sim.model$p),size=1)
}

# Rename before saving
ppmat.jp.abot  <- ppmat
constants.jp.abot  <- constants
d.jp.abot  <- d


save(constants.jp.abot,d.jp.abot,ppmat.jp.abot,file=here('results','ppcheck_jp_whole.RData'))
save(rhats,post.sample.core.jp.abot,file=here('results','post_jp_whole.RData'))
