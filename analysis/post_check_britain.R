# Load libraries, source code, data ----
library(here)
library(rcarbon)
source(here('src','utility.R'))

# Britian ----
load(here('results','post_gbwhole.RData'))
load(here('results','ppcheck_gbwhole.RData'))
load(here('results','post_gbregional.RData'))
load(here('results','ppcheck_gbregional.RData'))

# Extract relevant subsets ----
caldates.gbwhole  <- calibrate(d.gbwhole$cra,constants.gbwhole$cra_error)
obs.bin.gbwhole  <- d.gbwhole$y

caldates.gbregional  <- calibrate(d.gbregional$cra,constants.gbregional$cra_error)
obs.bin.gbregional  <- d.gbregional$y

caldates.engwales  <- caldates.gbregional[constants.gbregional$region==1]
caldates.scotmain  <- caldates.gbregional[constants.gbregional$region==2]
caldates.scotisles  <- caldates.gbregional[constants.gbregional$region==3]
obs.bin.engwales  <- obs.bin.gbregional[constants.gbregional$region==1]
obs.bin.scotmain  <- obs.bin.gbregional[constants.gbregional$region==2]
obs.bin.scotisles  <- obs.bin.gbregional[constants.gbregional$region==3]
ppmat.engwales  <- ppmat.gbregional[constants.gbregional$region==1,]
ppmat.scotmain  <- ppmat.gbregional[constants.gbregional$region==2,]
ppmat.scotisles  <- ppmat.gbregional[constants.gbregional$region==3,]

ppc.gbwhole  <- ppcheck(caldates.gbwhole,obs.bin.gbwhole,ppmat.gbregional,timeRange=c(8000,3000))
ppc.engwales  <- ppcheck(caldates.engwales,obs.bin.engwales,ppmat.engwales,timeRange=c(8000,3000))
ppc.scotmain  <- ppcheck(caldates.scotmain,obs.bin.scotmain,ppmat.scotmain,timeRange=c(8000,3000))
ppc.scotisles  <- ppcheck(caldates.scotisles,obs.bin.scotisles,ppmat.scotisles,timeRange=c(8000,3000))

save(ppc.gbwhole,file=here('results','ppc_gbwhole.RData'))
save(ppc.engwales,file=here('results','ppc_engwales.RData'))
save(ppc.scotmain,file=here('results','ppc_scotmain.RData'))
save(ppc.scotisles,file=here('results','ppc_scotisles.RData'))
