# Load Libraries & Results ----
library(here)
library(nimbleCarbon)
library(rcarbon)
source(here('src','utility.R'))


# Britain ----
load(here('results','post_gbwhole.RData'))
load(here('results','post_gbregional.RData'))
load(here('results','ppc_gbwhole.RData'))
load(here('results','ppc_engwales.RData'))
load(here('results','ppc_scotisles.RData'))
load(here('results','ppc_scotmain.RData'))

# Fitted model
par(mfrow=c(2,2))
plot.fitted(timeRange=c(8000,3000),r=post.sample.core.gbwhole[,'r'],m=post.sample.core.gbwhole[,'m'],mu_k=post.sample.core.gbwhole[,'mu_k'],calendar='BCAD',nsample=500)
title('Britain (whole)')

plot.fitted(timeRange=c(8000,3000),r=post.sample.core.gbregional[,'r[1]'],m=post.sample.core.gbregional[,'m[1]'],mu_k=post.sample.core.gbregional[,'mu_k[1]'],calendar='BCAD',nsample=500)
title('England and Wales')

plot.fitted(timeRange=c(8000,3000),r=post.sample.core.gbregional[,'r[2]'],m=post.sample.core.gbregional[,'m[2]'],mu_k=post.sample.core.gbregional[,'mu_k[2]'],calendar='BCAD',nsample=500)
title('Scotland (mainland)')

plot.fitted(timeRange=c(8000,3000),r=post.sample.core.gbregional[,'r[3]'],m=post.sample.core.gbregional[,'m[3]'],mu_k=post.sample.core.gbregional[,'mu_k[3]'],calendar='BCAD',nsample=500)
title('Scottish Isles')

# Posterior r
par(mfrow=c(2,2))
postHPDplot(post.sample.core.gbwhole[,'r'],xlab='r',ylab='Posterior Density',main='Britain (Whole)')
postHPDplot(post.sample.core.gbregional[,'r[1]'],xlab='r',ylab='Posterior Density',main='England and Wales')
postHPDplot(post.sample.core.gbregional[,'r[2]'],xlab='r',ylab='Posterior Density',main='Scotland (Mainland)')
postHPDplot(post.sample.core.gbregional[,'r[3]'],xlab='r',ylab='Posterior Density',main='Scottish Isles')

# Posterior m
par(mfrow=c(2,2))
postHPDplot(post.sample.core.gbwhole[,'m'],xlab='m',ylab='Posterior Density',main='Britain (Whole)')
postHPDplot(post.sample.core.gbregional[,'m[1]'],xlab='m',ylab='Posterior Density',main='England and Wales')
postHPDplot(post.sample.core.gbregional[,'m[2]'],xlab='m',ylab='Posterior Density',main='Scotland (Mainland)')
postHPDplot(post.sample.core.gbregional[,'m[3]'],xlab='m',ylab='Posterior Density',main='Scottish Isles')

# Predictive Checks
par(mfrow=c(2,2))
plotPcheck(ppc.gbwhole,calendar='BCAD')
title('Britain (whole)')
plotPcheck(ppc.engwales,calendar='BCAD')
title('England and Wales')
plotPcheck(ppc.scotmain,calendar='BCAD')
title('Scotland (mainland)')
plotPcheck(ppc.scotisles,calendar='BCAD')
title('Scottish Isles')


# Japan ----
load(here('results','post_jpwhole.RData'))
load(here('results','post_jpregional.RData'))
load(here('results','ppc_jpwhole.RData'))
load(here('results','ppc_central.RData'))
load(here('results','ppc_neast.RData'))
load(here('results','ppc_swest.RData'))

# Fitted model
par(mfrow=c(2,2))
plot.fitted(timeRange=c(4000,1700),r=post.sample.core.jpwhole[,'r'],m=post.sample.core.jpwhole[,'m'],mu_k=post.sample.core.jpwhole[,'mu_k'],calendar='BCAD',nsample=500)
title('Japan (whole)')

plot.fitted(timeRange=c(4000,1700),r=post.sample.core.jpregional[,'r[1]'],m=post.sample.core.jpregional[,'m[1]'],mu_k=post.sample.core.jpregional[,'mu_k[1]'],calendar='BCAD',nsample=500)
title('SW Japan')

plot.fitted(timeRange=c(4000,1700),r=post.sample.core.jpregional[,'r[2]'],m=post.sample.core.jpregional[,'m[2]'],mu_k=post.sample.core.jpregional[,'mu_k[2]'],calendar='BCAD',nsample=500)
title('C Japan')

plot.fitted(timeRange=c(4000,1700),r=post.sample.core.jpregional[,'r[3]'],m=post.sample.core.jpregional[,'m[3]'],mu_k=post.sample.core.jpregional[,'mu_k[3]'],calendar='BCAD',nsample=500)
title('NE Japan')

# Posterior r
par(mfrow=c(2,2))
postHPDplot(post.sample.core.jpwhole[,'r'],xlab='r',ylab='Posterior Density',main='Japan (whole)')
postHPDplot(post.sample.core.jpregional[,'r[1]'],xlab='r',ylab='Posterior Density',main='SW Japan')
postHPDplot(post.sample.core.jpregional[,'r[2]'],xlab='r',ylab='Posterior Density',main='Scotland (Mainland)')
postHPDplot(post.sample.core.jpregional[,'r[3]'],xlab='r',ylab='Posterior Density',main='NE Japan')

# Posterior m
par(mfrow=c(2,2))
postHPDplot(post.sample.core.jpwhole[,'m'],xlab='m',ylab='Posterior Density',main='Britain (Whole)')
postHPDplot(post.sample.core.jpregional[,'m[1]'],xlab='m',ylab='Posterior Density',main='SW Japan')
postHPDplot(post.sample.core.jpregional[,'m[2]'],xlab='m',ylab='Posterior Density',main='Scotland (Mainland)')
postHPDplot(post.sample.core.jpregional[,'m[3]'],xlab='m',ylab='Posterior Density',main='NE Japan')

# Predictive Checks
par(mfrow=c(2,2))
plotPcheck(ppc.jpwhole,calendar='BCAD')
title('Japan (whole)')
plotPcheck(ppc.swest,calendar='BCAD')
title('SW Japan')
plotPcheck(ppc.central,calendar='BCAD')
title('C Japan')
plotPcheck(ppc.neast,calendar='BCAD')
title('NE Japan')
