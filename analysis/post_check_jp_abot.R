# Load libraries, source code, data ----
library(here)
library(rcarbon)
source(here('src','utility.R'))

load(here('results','post_jp_abot.RData'))
load(here('results','ppcheck_jp_abot.RData'))

caldates.jp.abot  <- calibrate(d.jp.abot$cra,constants.jp.abot$cra_error)
obs.bin.jp.abot  <- d.jp.abot$y
fitted.mat  <- apply(ppmat.params.abot,1,function(x){sigmoid(4000:1700,r=x[1],m=x[2],k=x[3])})
fitted.mean.abot  <- apply(fitted.mat,1,mean)
fitted.90env.abot  <- apply(fitted.mat,1,quantile,prob=c(0.1,0.9))


ppc.jp.abot  <- ppcheck(caldates.jp.abot,obs.bin.jp.abot,ppmat.jp.abot,timeRange=c(4000,1700))

save(ppc.jp.abot,file=here('results','ppc_jp_abot.RData'))
