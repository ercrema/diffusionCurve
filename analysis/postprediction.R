# Load libraries, source code, data ----
library(here)
library(rcarbon)
source(here('src','utility.R'))
load(here('analysis','res_rice_3reg_core.RData'))
load(here('analysis','res_rice_3reg_ppcheck.RData'))

# Extract relevant subsets ----
caldates  <- calibrate(d$cra,constants$cra_error)
obs.bin  <- d$y

sw.caldates <- caldates[constants$region==1]
c.caldates <- caldates[constants$region==2]
ne.caldates <- caldates[constants$region==3]

sw.obs.bin  <- obs.bin[constants$region==1]
c.obs.bin  <- obs.bin[constants$region==2]
ne.obs.bin  <- obs.bin[constants$region==3]

sw.ppdata  <- ppmat[constants$region==1,]
c.ppdata  <- ppmat[constants$region==2,]
ne.ppdata  <- ppmat[constants$region==3,]

ppc.sw  <- ppcheck(sw.caldates,sw.obs.bin,sw.ppdata,timeRange=c(4000,1750))
ppc.c  <- ppcheck(c.caldates,c.obs.bin,c.ppdata,timeRange=c(4000,1750))
ppc.ne  <- ppcheck(ne.caldates,ne.obs.bin,ne.ppdata,timeRange=c(4000,1750))

save(ppc.sw,ppc.c,ppc.ne,file=here('analysis','posteriorchecks.RData'))
