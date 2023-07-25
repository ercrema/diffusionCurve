# Load libraries, source code, data ----
library(here)
library(rcarbon)
source(here('src','utility.R'))

# Japan ----
load(here('results','post_jpregional2.RData'))
load(here('results','ppcheck_jpregional2.RData'))

# Extract relevant subsets ----
caldates.jpregional2  <- calibrate(d.jpregional2$cra,constants.jpregional2$cra_error)
obs.bin.jpregional2  <- d.jpregional2$y

caldates.nkyushu  <- caldates.jpregional2[constants.jpregional2$region==1]
caldates.kinki  <- caldates.jpregional2[constants.jpregional2$region==2]
caldates.skanto  <- caldates.jpregional2[constants.jpregional2$region==3]
caldates.tohoku  <- caldates.jpregional2[constants.jpregional2$region==4]

obs.bin.nkyushu  <- obs.bin.jpregional2[constants.jpregional2$region==1]
obs.bin.kinki  <- obs.bin.jpregional2[constants.jpregional2$region==2]
obs.bin.skanto  <- obs.bin.jpregional2[constants.jpregional2$region==3]
obs.bin.tohoku  <- obs.bin.jpregional2[constants.jpregional2$region==4]

ppmat.nkyushu  <- ppmat.jpregional2[constants.jpregional2$region==1,]
ppmat.kinki  <- ppmat.jpregional2[constants.jpregional2$region==2,]
ppmat.skanto  <- ppmat.jpregional2[constants.jpregional2$region==3,]
ppmat.tohoku  <- ppmat.jpregional2[constants.jpregional2$region==4,]

ppc.nkyushu  <- ppcheck(caldates.nkyushu,obs.bin.nkyushu,ppmat.nkyushu,timeRange=c(4000,1700))
ppc.kinki  <- ppcheck(caldates.kinki,obs.bin.kinki,ppmat.kinki,timeRange=c(4000,1700))
ppc.skanto  <- ppcheck(caldates.skanto,obs.bin.skanto,ppmat.skanto,timeRange=c(4000,1700))
ppc.tohoku  <- ppcheck(caldates.tohoku,obs.bin.tohoku,ppmat.tohoku,timeRange=c(4000,1700))

save(ppc.nkyushu,file=here('results','ppc_nkyushu.RData'))
save(ppc.kinki,file=here('results','ppc_kinki.RData'))
save(ppc.skanto,file=here('results','ppc_skanto.RData'))
save(ppc.tohoku,file=here('results','ppc_tohoku.RData'))
