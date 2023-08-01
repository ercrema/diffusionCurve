# Load libraries, source code, data ----
library(here)
library(rcarbon)
source(here('src','utility.R'))

load(here('results','post_gb_abot.RData'))
load(here('results','ppcheck_gb_abot.RData'))

caldates.gb.abot  <- calibrate(d.gb.abot$cra,constants.gb.abot$cra_error)
obs.bin.gb.abot  <- d.gb.abot$y
fitted.mat  <- apply(ppmat.params.gb.abot,1,function(x){sigmoid(7000:3000,r=x[1],m=x[2],k=x[3])})
fitted.mean.gb.abot  <- apply(fitted.mat,1,mean)
fitted.90env.gb.abot  <- apply(fitted.mat,1,quantile,prob=c(0.1,0.9))
    
ppc.gb.abot  <- ppcheck(caldates.gb.abot,obs.bin.gb.abot,ppmat.gb.abot,timeRange=c(7000,3000))

save(ppc.gb.abot,fitted.mean.gb.abot,fitted.90env.gb.abot,file=here('results','ppc_gb_abot.RData'))
