# Load libraries, source code, data ----
library(here)
library(rcarbon)
source(here('src','utility.R'))

# Japan ----
load(here('results','post_jpwhole.RData'))
load(here('results','ppcheck_jpwhole.RData'))
load(here('results','post_jpregional.RData'))
load(here('results','ppcheck_jpregional.RData'))

# Extract relevant subsets ----
caldates.jpwhole  <- calibrate(d.jpwhole$cra,constants.jpwhole$cra_error)
obs.bin.jpwhole  <- d.jpwhole$y

caldates.jpregional  <- calibrate(d.jpregional$cra,constants.jpregional$cra_error)
obs.bin.jpregional  <- d.jpregional$y

caldates.swest  <- caldates.jpregional[constants.jpregional$region==1]
caldates.central  <- caldates.jpregional[constants.jpregional$region==2]
caldates.neast  <- caldates.jpregional[constants.jpregional$region==3]
obs.bin.swest  <- obs.bin.jpregional[constants.jpregional$region==1]
obs.bin.central  <- obs.bin.jpregional[constants.jpregional$region==2]
obs.bin.neast  <- obs.bin.jpregional[constants.jpregional$region==3]
ppmat.swest  <- ppmat.jpregional[constants.jpregional$region==1,]
ppmat.central  <- ppmat.jpregional[constants.jpregional$region==2,]
ppmat.neast  <- ppmat.jpregional[constants.jpregional$region==3,]

ppc.jpwhole  <- ppcheck(caldates.jpwhole,obs.bin.jpwhole,ppmat.jpwhole,timeRange=c(4000,1700))
ppc.swest  <- ppcheck(caldates.swest,obs.bin.swest,ppmat.swest,timeRange=c(4000,1700))
ppc.central  <- ppcheck(caldates.central,obs.bin.central,ppmat.central,timeRange=c(4000,1700))
ppc.neast  <- ppcheck(caldates.neast,obs.bin.neast,ppmat.neast,timeRange=c(4000,1700))

save(ppc.jpwhole,file=here('results','ppc_jpwhole.RData'))
save(ppc.swest,file=here('results','ppc_swest.RData'))
save(ppc.central,file=here('results','ppc_central.RData'))
save(ppc.neast,file=here('results','ppc_neast.RData'))
