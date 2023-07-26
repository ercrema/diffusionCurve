# Load libraries, source code, data ----
library(here)
library(rcarbon)
source(here('src','utility.R'))

load(here('results','post_gb_abot.RData'))
load(here('results','ppcheck_gb_abot.RData'))

caldates.gb.abot  <- calibrate(d.gb.abot$cra,constants.gb.abot$cra_error)
obs.bin.gb.abot  <- d.gb.abot$y
ppc.gb.abot  <- ppcheck(caldates.gb.abot,obs.bin.gb.abot,ppmat.gb.abot,timeRange=c(7000,3000))

save(ppc.gb.abot,file=here('results','ppc_gb_abot.RData'))
