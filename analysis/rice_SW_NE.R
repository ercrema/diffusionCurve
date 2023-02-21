# Load libraries and Data ----
library(here)
library(parallel)
library(nimbleCarbon)
library(rcarbon)
data(intcal20)
load(here('data','ricedata.RData'))

# Data Preparation ----
caldates  <- calibrate(ricedata$C14Age,ricedata$C14Error)
ii  <- which.CalDates(caldates,BP<=4000&BP>1700,p=0.5)
ricedata <- ricedata[ii,]
caldates  <- caldates[ii]
ricedata$site.id  <- as.integer(factor(ricedata$SiteNameJp))
ricedata$Region2  <- ifelse(ricedata$Region %in% c('Chubu','Kanto','Tohoku'),1,2)


# Define Constants ----
constants  <- list()
constants$N  <- nrow(ricedata)
constants$NSites  <- length(unique(ricedata$site.id))
constants$siteID  <- ricedata$site.id
constants$region  <- ricedata$Region2
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma
constants$cra_error  <- ricedata$C14Error

# Define Data ----
d  <- list()
d$y  <- ricedata$rice_nuts 
d$cra  <- ricedata$C14Age

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
			p[i]  <- k[siteID[i]]/(1+exp(r[region[i]]*(theta[i]-m[region[i]]))); #sigmoidal model


			#Calibration of observed cra
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sigma[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sigma[i]);

			# Prior for theta (flat)
			theta[i] ~ dunif(1000,10000)
		}


		r[1] ~ dexp(10) # prior adoption rate
		r[2] ~ dexp(10) # prior adoption rate
		m[1] ~ T(dnorm(mean=2500,sd=500),1000,50000) #prior mid-point
		m[2] ~ T(dnorm(mean=2500,sd=500),1000,50000) #prior mid-point

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
	inits$r  <- c(0.0001,0.0001)
	inits$m  <- c(3000,3000)
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
out  <- parLapply(cl=cl,X=seeds,fun=runFun,d=d,constants=constants,theta=theta.init,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(cl)

# Diagnostic and Posterior Processing ----
post.sample  <- coda::mcmc.list(out)
diagnostic  <- coda::gelman.diag(post.sample)

