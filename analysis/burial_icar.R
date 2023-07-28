library(here)
library(nimble)
library(rcarbon)
library(nimbleCarbon)
library(parallel)
source(here('src','utility.R'))
load(here('data','burialdata.RData'))

# Define Constants and Data
icar.struct  <- function(x)
{
	num  <- c(1,rep(2,x-2),1)
	adj  <- numeric()
	for (i in 1:x)
	{
		if(i==1){adj = 2}
		if(i==x){adj = c(adj,x-1)}
		if(i>1&i<x){adj = c(adj,c(i-1,i+1))}
	}
	weight <- rep(1,length(adj))
	return(list(adj=adj,weight=weight,num=num))
}

data(intcal20)
constants  <- list()
constants$n  <- nrow(burial)
constants$res  <- 100 #100 year resolution
constants$a  <- 5500 #start date
constants$b  <- 2000 #end 
constants$n.tblocks  <- 1 + (constants$a - constants$b)/constants$res 
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma
constants$weights <- icar.struct(constants$n.tblocks)$weight
constants$num <- icar.struct(constants$n.tblocks)$num
constants$adj <- icar.struct(constants$n.tblocks)$adj
constants$L  <- length(constants$adj)

# Setup Data ----
d <- list(y=burial$Cremation,cra=burial$CRA,cra_error=burial$Error)
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
post.sample  <- coda::mcmc.list(out)
rhats  <- coda::gelman.diag(post.sample)
# which(rhats[[1]][,1]>1.01)
# rhats[[1]][which(rhats[[1]][,1]>1.01),1]
# lower convergence only for theta parameters

# Store output ----
post.sample.combined.icar  <- do.call(rbind.data.frame,post.sample)
agr  <- agreementIndex(d$cra,d$cra_error,theta=post.sample.combined.icar[,grep('theta',colnames(post.sample.combined.icar))])
min(agr$agreement)
# theta show good agreement index, with smallest value equal to 72, indicating that the model per se is not problematic
post.sample.combined.icar  <- post.sample.combined.icar[,grep('pseq',colnames(post.sample.combined.icar))]
constants.icar  <- constants
rhats.icar  <- rhats
save(rhats.icar,post.sample.combined.icar,constants.icar,file=here('results','post_icar_burial.RData'))

