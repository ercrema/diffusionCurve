library(nimbleCarbon)
library(parallel)
library(rcarbon)
library(here)
load(here('sim','simdata','simdata1a.RData'))

# Remove true.parameters from constants
constants$mu_k  <- NULL
constants$phi  <- NULL

# Inits
caldates  <- calibrate(d$cra,constants$cra_error)
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
rhats.sim1a  <- coda::gelman.diag(post.sample)
which(rhats.sim1a[[1]][,1]>1.01) #only thetas
rhats.sim1a[[1]][1:4,]

# Store output ----
post.sample.combined  <- do.call(rbind.data.frame,post.sample)
post.sample.theta  <- post.sample.combined[,grep('theta',colnames(post.sample.combined))]
post.sample.core.sim1a  <- post.sample.combined[,!grepl('theta',colnames(post.sample.combined))]

save(rhats.sim1a,post.sample.core.sim1a,file=here('sim','results','post_sim1a.RData'))
