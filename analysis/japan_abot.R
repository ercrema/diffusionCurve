# Load libraries and Data ----
library(here)
library(parallel)
library(nimbleCarbon)
library(rcarbon)
data(intcal20)
load(here('data','jpdata.RData'))

# Data Preparation ----
caldates  <- calibrate(jpdata$C14Age,jpdata$C14Error)
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
			k[j] ~ dbeta(beta0,beta1)
		}

		mu_k ~ dbeta(2,2) 
		phi ~ dgamma(5,0.1)
		beta0  <- mu_k * (phi) + 1
		beta1  <- (1 - mu_k) * (phi) + 1
	})

	#Define inits
	inits  <- list()
	inits$r  <- 0.0001
	inits$m  <- 3000
	inits$mu_k  <- 0.5
	inits$phi  <- 10
	inits$k  <- rbeta(constants$NSites,(inits$mu_k*(inits$phi) +1),((1-inits$mu_k)*(inits$phi)+1))
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
which(rhats.jp.abot[[1]][,1]>1.01) #only thetas

# Store output ----
post.sample.combined  <- do.call(rbind.data.frame,post.sample)
post.sample.theta  <- post.sample.combined[,grep('theta',colnames(post.sample.combined))]
post.sample.core.jp.abot  <- post.sample.combined[,!grepl('theta',colnames(post.sample.combined))]

# Agreement index for theta ----
agree.jp.abot <- agreementIndex(CRA=d$cra,CRAError=constants$cra_error,theta=post.sample.theta)
# min(agree.jp.abot$agreement) 90

# Posterior Predictive Checks ----
nsim  <- 1000 #Number of posterior simulations
ppmat  <- matrix(NA,nrow=constants$N,ncol=nsim) #Matrix storing predictions
ppmat.params  <- matrix(NA,ncol=4,nrow=nsim) |> as.data.frame()
colnames(ppmat.params)  <- c('r','m','mu_k','phi')
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
		k[j] ~ dbeta(beta0,beta1)
	}

	mu_k ~ dbeta(2,2)
	phi ~ dgamma(5,0.1)
	beta0  <- mu_k * (phi) + 1
	beta1  <- (1 - mu_k) * (phi) + 1
})

s.index  <- sample((niter-nburnin)*nchains/thin,size=nsim)
ppmat.params$phi  <- as.numeric(post.sample.combined[s.index,'phi'])
ppmat.params$mu_k  <- as.numeric(post.sample.combined[s.index,'mu_k'])
ppmat.params$r  <- as.numeric(post.sample.combined[s.index,'r'])
ppmat.params$m  <- as.numeric(post.sample.combined[s.index,'m'])

sim.model  <- nimbleModel(adoptionSimModel,constants=constants,data=d)

pb <- txtProgressBar(min = 0, max = nsim, style = 3, width = 50, char = "=")

for (i in 1:nsim)
{
    setTxtProgressBar(pb, i)
    ii  <- s.index[i]
    sim.model$phi  <- as.numeric(post.sample.combined[ii,'phi'])
    sim.model$mu_k  <- as.numeric(post.sample.combined[ii,'mu_k'])
    sim.model$r  <- as.numeric(post.sample.combined[ii,'r'])
    sim.model$m  <- as.numeric(post.sample.combined[ii,'m'])
    sim.model$theta  <- as.numeric(post.sample.combined[ii,grep('theta\\[',colnames(post.sample.combined))])
    sim.model$calculate('beta0')
    sim.model$calculate('beta1')
    sim.model$simulate('k')
    sim.model$simulate('p')
    ppmat[,i]  <- rbinom(constants$N,prob=unlist(sim.model$p),size=1)
}

# Rename before saving
ppmat.jp.abot  <- ppmat
ppmat.params.jp.abot  <- ppmat.params
constants.jp.abot  <- constants
d.jp.abot  <- d


save(constants.jp.abot,d.jp.abot,ppmat.jp.abot,ppmat.params.jp.abot,file=here('results','ppcheck_jp_abot.RData'))
save(agree.jp.abot,rhats.jp.abot,post.sample.core.jp.abot,file=here('results','post_jp_abot.RData'))
