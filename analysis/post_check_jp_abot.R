# Load libraries, source code, data ----
library(here)
library(rcarbon)
source(here('src','utility.R'))

load(here('results','post_jp_abot.RData'))
load(here('results','ppcheck_jp_abot.RData'))

caldates.jp.abot  <- calibrate(d.jp.abot$cra,constants.jp.abot$cra_error)
obs.bin.jp.abot  <- d.jp.abot$y

ppc.jp.abot  <- ppcheck(caldates.jp.abot,obs.bin.jp.abot,ppmat.jp.abot,timeRange=c(4000,1700))

save(ppc.jp.abot,file=here('results','ppc_jp_abot.RData'))
