# Load Libraries & Results ----
library(here)
library(ggplot2)
library(tidybayes)
library(tidyr)
library(gridExtra)
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
pdf(here('figures','fitted_britain.pdf'),width=6,height=10)
par(mfrow=c(3,1))
# plot.fitted(timeRange=c(8000,3000),r=post.sample.core.gbwhole[,'r'],m=post.sample.core.gbwhole[,'m'],mu_k=post.sample.core.gbwhole[,'mu_k'],calendar='BCAD',nsample=500)
# title('Britain (whole)')
plot.fitted(timeRange=c(8000,3000),r=post.sample.core.gbregional[,'r[1]'],m=post.sample.core.gbregional[,'m[1]'],mu_k=post.sample.core.gbregional[,'mu_k[1]'],calendar='BCAD',nsample=500)
title('England and Wales')

plot.fitted(timeRange=c(8000,3000),r=post.sample.core.gbregional[,'r[2]'],m=post.sample.core.gbregional[,'m[2]'],mu_k=post.sample.core.gbregional[,'mu_k[2]'],calendar='BCAD',nsample=500)
title('Scotland (mainland)')

plot.fitted(timeRange=c(8000,3000),r=post.sample.core.gbregional[,'r[3]'],m=post.sample.core.gbregional[,'m[3]'],mu_k=post.sample.core.gbregional[,'mu_k[3]'],calendar='BCAD',nsample=500)
title('Scottish Isles')
dev.off()

# Posterior r
post.r.gbregional  <- post.sample.core.gbregional[,grep('r\\[',colnames(post.sample.core.gbregional))]
colnames(post.r.gbregional)  <- c('England and Wales','Scotland (mainland)','Scottish Isles')
post.r.gbregional  <- gather(post.r.gbregional)

# pdf(here('figures','post_r_britain.pdf'),width=7,height=6)
r.britain <- ggplot(aes(x = key, y = value),data=post.r.gbregional) + 
	     stat_interval(aes(x = key), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	     scale_color_brewer() +
	     xlab('')   +
	     ylab('Adoption Rate') +
	     ggtitle('Estimated Adoption Rate (Oat+Barley/Hazelnut) in Britain (8000-3000 BP)') +
	     theme(axis.text.y = element_text(angle = 90))
# dev.off()

# Posterior mu_k
post.k.gbregional  <- post.sample.core.gbregional[,grep('mu_k\\[',colnames(post.sample.core.gbregional))]
colnames(post.k.gbregional)  <- c('England and Wales','Scotland (mainland)','Scottish Isles')
post.k.gbregional  <- gather(post.k.gbregional)
post.k.gbregional$value  <- logistic(post.k.gbregional$value)

# pdf(here('figures','post_k_britain.pdf'),width=7,height=6)
k.britain  <- ggplot(aes(x = key, y = value),data=post.k.gbregional) + 
	      stat_interval(aes(x = key), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	      scale_color_brewer() +
	      xlab('')   +
	      ylab('Adoption Cap') +
	      ggtitle('Estimated Adoption Cap (Oat+Barley/Hazelnut) in Britain (8000-3000 BP)') +
	      theme(axis.text.y = element_text(angle = 90))
# dev.off()     

# Posterior m
post.m.gbregional  <- post.sample.core.gbregional[,grep('m\\[',colnames(post.sample.core.gbregional))]
colnames(post.m.gbregional)  <- c('England and Wales','Scotland (mainland)','Scottish Isles')
post.m.gbregional  <- gather(post.m.gbregional)
post.m.gbregional$value = -post.m.gbregional$value #negative as tidybayes struggles in reversing x-axis

# pdf(here('figures','post_m_britain.pdf'),width=7,height=6)
m.britain <- ggplot(aes(y = value, x = key),data=post.m.gbregional) + 
	     stat_interval(aes(y = value), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	     scale_color_brewer() +
	     xlab('BP')   +
	     ylab('') +
	     ggtitle('Estimated Mid-point (Oat+Barley/Hazelnut) in Britain (8000-3000 BP)') +
# 	     theme(axis.text.x = element_text(angle = 90)) + 
	     scale_y_continuous(labels = abs)
# dev.off()

pdf(here('figures','post_rmk_britain.pdf'),width=6.5,height=10)
grid.arrange(r.britain,k.britain,m.britain,ncol=1)
dev.off()


# Predictive Checks
pdf(here('figures','pcheck_britain.pdf'),width=6,height=10)
par(mfrow=c(3,1))
# plotPcheck(ppc.gbwhole,calendar='BCAD')
# title('Britain (whole)')
plotPcheck(ppc.engwales,calendar='BCAD')
title('England and Wales')
plotPcheck(ppc.scotmain,calendar='BCAD')
title('Scotland (mainland)')
plotPcheck(ppc.scotisles,calendar='BCAD')
title('Scottish Isles')
dev.off()

# Japan ----
load(here('results','post_jpwhole.RData'))
load(here('results','post_jpregional.RData'))
load(here('results','ppc_jpwhole.RData'))
load(here('results','ppc_central.RData'))
load(here('results','ppc_neast.RData'))
load(here('results','ppc_swest.RData'))

# Fitted model
pdf(here('figures','fitted_japan.pdf'),width=6,height=10)
par(mfrow=c(3,1))
# plot.fitted(timeRange=c(4000,1700),r=post.sample.core.jpwhole[,'r'],m=post.sample.core.jpwhole[,'m'],mu_k=post.sample.core.jpwhole[,'mu_k'],calendar='BCAD',nsample=500)
# title('Japan (whole)')

plot.fitted(timeRange=c(4000,1700),r=post.sample.core.jpregional[,'r[1]'],m=post.sample.core.jpregional[,'m[1]'],mu_k=post.sample.core.jpregional[,'mu_k[1]'],calendar='BCAD',nsample=500)
title('SW Japan')

plot.fitted(timeRange=c(4000,1700),r=post.sample.core.jpregional[,'r[2]'],m=post.sample.core.jpregional[,'m[2]'],mu_k=post.sample.core.jpregional[,'mu_k[2]'],calendar='BCAD',nsample=500)
title('C Japan')

plot.fitted(timeRange=c(4000,1700),r=post.sample.core.jpregional[,'r[3]'],m=post.sample.core.jpregional[,'m[3]'],mu_k=post.sample.core.jpregional[,'mu_k[3]'],calendar='BCAD',nsample=500)
title('NE Japan')
dev.off()


# Posterior r
post.r.jpregional  <- post.sample.core.jpregional[,grep('r\\[',colnames(post.sample.core.jpregional))]
colnames(post.r.jpregional)  <- c('South-West','Central','North-East')
post.r.jpregional  <- gather(post.r.jpregional)
post.r.jpregional$key  <- factor(post.r.jpregional$key,c('South-West','Central','North-East'))

# pdf(here('figures','post_r_japan.pdf'),width=7,height=6)
r.japan <- ggplot(aes(x = key, y = value),data=post.r.jpregional) + 
	   stat_interval(aes(x = key), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	   scale_color_brewer() +
	   xlab('')   +
	   ylab('Adoption Mid-point') +
	   ggtitle('Estimated Adoption Mid-point (Rice/Nuts) in Japan (2000BC-300CE)') +
	   theme(axis.text.y = element_text(angle = 90))
# dev.off()


# Posterior k
post.k.jpregional  <- post.sample.core.jpregional[,grep('mu_k\\[',colnames(post.sample.core.jpregional))]
colnames(post.k.jpregional)  <- c('South-West','Central','North-East')
post.k.jpregional  <- gather(post.k.jpregional)
post.k.jpregional$key  <- factor(post.k.jpregional$key,c('South-West','Central','North-East'))
post.k.jpregional$value  <- logistic(post.k.jpregional$value)


# pdf(here('figures','post_k_japan.pdf'),width=7,height=6)
k.japan <- ggplot(aes(x = key, y = value),data=post.k.jpregional) + 
	   stat_interval(aes(x = key), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	   scale_color_brewer() +
	   xlab('')   +
	   ylab('Adoption Cap') +
	   ggtitle('Estimated Adoption Cap (Rice/Nuts) in Japan (2000BC-300CE)') +
	   theme(axis.text.y = element_text(angle = 90))
# dev.off()


# pdf(here('figures','post_r_japan.pdf'),width=8,height=8)
# par(mfrow=c(2,2))
# postHPDplot(post.sample.core.jpwhole[,'r'],xlab='r',ylab='Posterior Density',main='Japan (whole)',xlim=c(0,0.3))
# postHPDplot(post.sample.core.jpregional[,'r[1]'],xlab='r',ylab='Posterior Density',main='SW Japan',xlim=c(0,0.3))
# postHPDplot(post.sample.core.jpregional[,'r[2]'],xlab='r',ylab='Posterior Density',main='Scotland (Mainland)',xlim=c(0,0.3))
# postHPDplot(post.sample.core.jpregional[,'r[3]'],xlab='r',ylab='Posterior Density',main='NE Japan',xlim=c(0,0.3))
# dev.off()
# 
# Posterior m
post.m.jpregional  <- post.sample.core.jpregional[,grep('m\\[',colnames(post.sample.core.jpregional))]
colnames(post.m.jpregional)  <- c('South-West','Central','North-East')
post.m.jpregional  <- gather(post.m.jpregional)
post.m.jpregional$key  <- factor(post.m.jpregional$key,c('South-West','Central','North-East'))
post.m.jpregional$value  <- BPtoBCAD(post.m.jpregional$value)

# pdf(here('figures','post_m_japan.pdf'),width=7,height=6)
m.japan  <- ggplot(aes(y = value, x = key),data=post.m.jpregional) + 
	    stat_interval(aes(x = key), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	    scale_color_brewer() +
	    xlab('BCE/CE')   +
	    ylab('') +
	    ggtitle('Estimated Adoption Midpoint (Rice/Nuts) in Japan (2000BC-300CE)') +
# 	    theme(axis.text.y = element_text(angle = 90)) +
	    scale_y_continuous(labels = abs)
# dev.off()

pdf(here('figures','post_rmk_japan.pdf'),width=6.5,height=10)
grid.arrange(r.japan,k.japan,m.japan,ncol=1)
dev.off()

# Predictive Checks
pdf(here('figures','pcheck_japan.pdf'),width=6,height=10)
par(mfrow=c(3,1))
# plotPcheck(ppc.jpwhole,calendar='BCAD')
# title('Japan (whole)')
plotPcheck(ppc.swest,calendar='BCAD')
title('SW Japan')
plotPcheck(ppc.central,calendar='BCAD')
title('C Japan')
plotPcheck(ppc.neast,calendar='BCAD')
title('NE Japan')
dev.off()

# Burial ----
load(here('results','post_icar_burial.RData'))


tblocks <- seq(constants.icar$a,constants.icar$b,by=-constants.icar$res)


pdf(here('figures','icar_cremation.pdf'),width=6,height=5)
plot(tblocks,apply(post.sample.combined.icar[,1:constants.icar$n.tblocks],2,mean),pch=20,xlim=c(constants.icar$a,constants.icar$b),ylim=c(0,1),type='n',xlab='Cal BP',ylab='Estimated Proportion of Cremations')
for (i in 1:constants.icar$n.tblocks)
{
	rect(xleft=tblocks[i]+45,xright=tblocks[i]-45,ybottom=quantile(post.sample.combined.icar[,i],0.025),ytop=quantile(post.sample.combined.icar[,i],0.975),border=NA,col='lightblue')
	rect(xleft=tblocks[i]+45,xright=tblocks[i]-45,ybottom=quantile(post.sample.combined.icar[,i],0.25),ytop=quantile(post.sample.combined.icar[,i],0.75),border=NA,col='steelblue')

}
lines(tblocks,apply(post.sample.combined.icar[,1:constants.icar$n.tblocks],2,mean),pch=20,type='b')
legend(x=3900,y=0.2,legend=c('95% HPD','50% HPD'),fill=c('lightblue','steelblue'),cex=0.9)
dev.off()
