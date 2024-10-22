library(here)
library(nimble)
library(rcarbon)
library(nimbleCarbon)
library(parallel)
source(here('src','utility.R'))
load(here('sim','simdata','simdata2.RData'))
theta.init  <- medCal(calibrate(d$cra,d$cra_error))
theta.init  <- ifelse(theta.init>constants$a-1,constants$a,theta.init)
theta.init  <- ifelse(theta.init<constants$b+2,constants$b,theta.init)


# Core runscript ----
runFun  <- function(seed,d,constants,theta,nburnin,niter,thin)
{
	#Load libraries
	library(nimbleCarbon)

	#Define model
	adoptionICAR  <- nimbleCode({
		for (i in 1:n)
		{
			y[i] ~ dbern(p[i]) #probability of event y
			p[i]  <- pseq[theta.index[i]]


			#Calibration of observed cra
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sigma[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sigma[i]);

			# Prior for theta (flat)
			theta[i] ~ dunif(b,a)

			# Round theta and assign to timeblock
			theta.rnd[i]  <- round(theta[i]/res)*res
			theta.index[i]  <- 1 + (a - theta.rnd[i])/res
		}

		for (j in 1:n.tblocks)
		{
			pseq[j] <- 1 /(1 + exp(-lpseq[j]))
		}

		lpseq[1:n.tblocks] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.tblocks], tau, zero_mean = 0)
		tau <- 1/sigmap^2
		sigmap  ~ dexp(1)

	})

	# Setup Inits
	set.seed(seed)
	inits  <- list()
	inits$theta  <- theta
	inits$sigma  <- rexp(1)
	inits$lpseq  <- rnorm(constants$n.tblocks,0,0.5)

	# Setup MCMC
	model  <- nimbleModel(adoptionICAR,constants=constants,inits=inits,data=d)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model)
	conf$addMonitors('pseq')
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)

	#run model
	results <- runMCMC(cMCMC, nchain=1,niter = niter, thin=thin, nburnin = nburnin,samplesAsCodaMCMC = T,setSeed=seed) 
}

# Parallelisation Setup ----
niter  <- 200000
nburnin <- 100000
thin  <- 4
nchains  <- 4
cl  <- makeCluster(nchains)
seeds  <- c(12,34,56,78)[1:nchains]

# Run MCMC ----
out  <- parLapply(cl=cl,X=seeds,fun=runFun,d=d,constants=constants,theta=theta.init,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(cl)

# Diagnostic and Posterior Processing ----
post.sample.sim2  <- coda::mcmc.list(out)
rhats.sim2  <- coda::gelman.diag(post.sample.sim2)
# rhats.sim3[[1]][1:36,1]
# Store output ----
post.sample.combined.icar.sim2  <- do.call(rbind.data.frame,post.sample.sim2)
post.sample.combined.icar.sim2  <- post.sample.combined.icar.sim2[,grep('pseq',colnames(post.sample.combined.icar.sim2))]
constants.icar.sim2  <- constants
save(post.sample.combined.icar.sim2,constants.icar.sim2,file=here('sim','results','post_icar_sim2.RData'))

