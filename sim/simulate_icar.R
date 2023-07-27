library(here)
library(nimble)
library(rcarbon)
library(nimbleCarbon)
source(here('src','utility.R'))
seed  <- 123

# Simulate Data ----
set.seed(seed)
timeRange  <- c(5500,2000)
BP  <- c(5000,4750,4500,4250,4000,3750,3500,3350,3300,2000)
logP  <- c(-0.5,-2,-1,0,-1,-0.8,-0.6,-0.2,0.3,0) 
fit  <- lm(logP~poly(BP,7))
Pseq  <- logistic(predict(fit,data.frame(BP=timeRange[1]:timeRange[2])))
# plot(timeRange[1]:timeRange[2],Pseq,xlim=timeRange,type='l')
n  <- 2536 #number of samples
caldates  <- round(runif(n,timeRange[2],timeRange[1]))
p  <- logistic(predict(fit,data.frame(BP=caldates)))
sim.d  <- data.frame(x=rbinom(n,1,prob=p),p=p,c14age=uncalibrate(caldates)[,4],calage=caldates)

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
constants$n  <- n
constants$res  <- 100 #100 year resolution
constants$a  <- 5500 #start date
constants$b  <- 2000 #end date
constants$n.tblocks  <- 1 + (constants$a - constants$b)/constants$res 
constants$calBP  <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err  <- intcal20$C14Age.sigma
constants$weights <- icar.struct(constants$n.tblocks)$weight
constants$num <- icar.struct(constants$n.tblocks)$num
constants$adj <- icar.struct(constants$n.tblocks)$adj
constants$L  <- length(constants$adj)

d <- list(y=sim.d$x,cra=sim.d$c14age,cra_error=rep(20,n))

save(Pseq,timeRange,d,constants,file=here('sim','simdata','simdata3.RData'))




# theta.init  <- medCal(calibrate(d$cra,d$cra_error))
# theta.init  <- ifelse(theta.init>constants$a-1,constants$a,theta.init)
# theta.init  <- ifelse(theta.init<constants$b+2,constants$b,theta.init)
# # Setup Inits ----
# inits  <- list()
# inits$theta  <- theta.init
# inits$sigma  <- 1
# inits$lpseq  <- rep(0,constants$n.tblocks)
# 
# 
# 
# # Define model ----
# adoptionICAR  <- nimbleCode({
# 	for (i in 1:n)
# 	{
# 		y[i] ~ dbern(p[i]) #probability of event y
# 		p[i]  <- pseq[theta.index[i]]
# 
# 
# 		#Calibration of observed cra
# 		mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
# 		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
# 		sigma[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
# 		cra[i] ~ dnorm(mean=mu[i],sd=sigma[i]);
# 
# 		# Prior for theta (flat)
# 		theta[i] ~ dunif(b,a)
# 
# 		# Round theta and assign to timeblock
# 		theta.rnd[i]  <- round(theta[i]/res)*res
# 		theta.index[i]  <- 1 + (a - theta.rnd[i])/res
# 	}
# 
# 	for (j in 1:n.tblocks)
# 	{
# 		pseq[j] <- 1 /(1 + exp(-lpseq[j]))
# 	}
# 
# 	lpseq[1:n.tblocks] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.tblocks], tau, zero_mean = 0)
# 	tau <- 1/sigmap^2
# 	sigmap  ~ dexp(1)
# 
# })
# 
# # Setup MCMC ----
# model  <- nimbleModel(adoptionICAR,constants=constants,inits=inits,data=d)
# cModel <- compileNimble(model)
# conf <- configureMCMC(model)
# conf$addMonitors('pseq')
# MCMC <- buildMCMC(conf)
# cMCMC <- compileNimble(MCMC)
# results <- runMCMC(cMCMC, nchain=2,niter = 10000, thin=1, nburnin =2000,samplesAsCodaMCMC = T) 
# res  <- do.call(rbind.data.frame,results)
# 
# tblocks <- seq(constants$a,constants$b,by=-constants$res)
# plot(tblocks,apply(res[,1:30],2,mean),pch=20,xlim=c(constants$a,constants$b),ylim=c(0,1))
# points(tblocks,apply(res[,1:30],2,quantile,0.025),pch="+",xlim=c(constants$a,constants$b),ylim=c(0,1))
# points(tblocks,apply(res[,1:30],2,quantile,0.975),pch="+",xlim=c(constants$a,constants$b),ylim=c(0,1))
# lines(constants$a:constants$b,logistic(predict(fit,data.frame(BP=constants$a:constants$b))))
# 
