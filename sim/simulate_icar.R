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

save(Pseq,timeRange,d,constants,file=here('sim','simdata','simdata2.RData'))
