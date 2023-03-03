# Load Libraries & Results ----
library(here)
library(nimbleCarbon)
library(rcarbon)
source(here('src','utility.R'))

load(here('analysis','posteriorchecks.RData'))
load(here('analysis','res_rice_3reg_core.RData'))

# Compare Overall Fit
pdf(width=10,height=7,file=here('analysis','fitted_model.pdf'))
par(mfrow=c(1,3))
plot.fitted(timeRange=c(4000,1750),r=post.sample.core[,'r[1]'],m=post.sample.core[,'m[1]'],mu_k=post.sample.core[,'mu_k[1]'],calendar='BCAD',nsample=500)
title('SW Japan')
plot.fitted(timeRange=c(4000,1750),r=post.sample.core[,'r[2]'],m=post.sample.core[,'m[2]'],mu_k=post.sample.core[,'mu_k[2]'],calendar='BCAD',nsample=500)
title('C Japan')
plot.fitted(timeRange=c(4000,1750),r=post.sample.core[,'r[3]'],m=post.sample.core[,'m[3]'],mu_k=post.sample.core[,'mu_k[3]'],calendar='BCAD',nsample=500)
title('NE Japan')
dev.off()

pdf(width=10,height=7,file=here('analysis','ppcheck.pdf'))
par(mfrow=c(1,3))
plotPcheck(ppc.sw ,calendar='BCAD')
title('SW Japan')
plotPcheck(ppc.c ,calendar='BCAD')
title('C Japan')
plotPcheck(ppc.ne ,calendar='BCAD')
title('NE Japan')
dev.off()



par(mfrow=c(1,3))
postHPDplot((post.sample.core[,'r[1]']*120))
title('SW Japan')
postHPDplot((post.sample.core[,'r[2]']*104))
title('C Japan')
postHPDplot((post.sample.core[,'r[3]']*66))
title('NE Japan')


par(mfrow=c(1,3))
postHPDplot(BPtoBCAD(post.sample.core[,'m[1]']))
title('SW Japan')
postHPDplot(BPtoBCAD(post.sample.core[,'m[2]']))
title('C Japan')
postHPDplot(BPtoBCAD(post.sample.core[,'m[3]']))
title('NE Japan')
