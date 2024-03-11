library(here)
library(nimbleCarbon)
library(latex2exp)
library(RColorBrewer)
library(emdbook)
library(coda)

# Figure S1 ----
load(here('sim','simdata','simdata1a.RData'))
load(here('sim','results','post_sim1a.RData'))

pdf(file=here('figures_and_tables','figureS1.pdf'),width=7,height=7)
layout(matrix(c(1,2,3,4,5,10,14,15,16,6,20,11,17,18,7,20,20,12,19,8,20,20,20,13,9),nrow=5,ncol=5),width=c(0.1,1,1,1,1),height=c(1,1,1,1,0.1))

label.plot <- function(label)
{
	par(mar=c(0,0,0,0))
	plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
	text(label,x=0.5,y=0.5,las=2,adj=0.5,cex=1.5)
}

label.plot('r')
label.plot('m')
label.plot(TeX(r'($\mu$)'))
label.plot(TeX(r'($\varphi$)'))

label.plot('')
label.plot('r')
label.plot('m')
label.plot(TeX(r'($\mu$)'))
label.plot(TeX(r'($\varphi$)'))

par(mar=c(2,2,2,2))
postHPDplot(post.sample.core.sim1a$r,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
abline(v=true.param.1a$r,lty=2,col=2)
postHPDplot(post.sample.core.sim1a$m,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
abline(v=true.param.1a$m,lty=2,col=2)
postHPDplot(post.sample.core.sim1a$mu_k,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
abline(v=true.param.1a$mu_k,lty=2,col=2)
postHPDplot(post.sample.core.sim1a$phi,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
abline(v=true.param.1a$phi,lty=2,col=2)


post.sample.core.sim1a  <- post.sample.core.sim1a[,c('r','m','mu_k','phi')]
true.param.1a <- true.param.1a[c('r','m','mu_k','phi')]


for (i in 1:4)
{
	for (j in 1:4)
	{
		if (j>i){
			plot(post.sample.core.sim1a[,i],post.sample.core.sim1a[,j],pch=1,col=rgb(0.67,0.84,0.9,0.1))
			HPDregionplot(mcmc(post.sample.core.sim1a[,c(i,j)]),add=TRUE,prob=c(0.5,0.95),n=100,lty=c(1,2))
			points(true.param.1a[i],true.param.1a[j],pch="+",col='red',cex=1.5)
		}
	}
}
dev.off()




# Figure S2 ----
load(here('sim','simdata','simdata1b.RData'))
load(here('sim','results','post_sim1b.RData'))

pdf(file=here('figures_and_tables','figureS2.pdf'),width=7,height=7)
layout(matrix(c(1,2,3,4,5,10,14,15,16,6,20,11,17,18,7,20,20,12,19,8,20,20,20,13,9),nrow=5,ncol=5),width=c(0.1,1,1,1,1),height=c(1,1,1,1,0.1))

label.plot <- function(label)
{
	par(mar=c(0,0,0,0))
	plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
	text(label,x=0.5,y=0.5,las=2,adj=0.5,cex=1.5)
}

label.plot('r')
label.plot('m')
label.plot(TeX(r'($\mu$)'))
label.plot(TeX(r'($\varphi$)'))

label.plot('')
label.plot('r')
label.plot('m')
label.plot(TeX(r'($\mu$)'))
label.plot(TeX(r'($\varphi$)'))

par(mar=c(2,2,2,2))
postHPDplot(post.sample.core.sim1b$r,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
abline(v=true.param.1b$r,lty=2,col=2)
postHPDplot(post.sample.core.sim1b$m,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
abline(v=true.param.1b$m,lty=2,col=2)
postHPDplot(post.sample.core.sim1b$mu_k,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
abline(v=true.param.1b$mu_k,lty=2,col=2)
postHPDplot(post.sample.core.sim1b$phi,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
abline(v=true.param.1b$phi,lty=2,col=2)


post.sample.core.sim1b  <- post.sample.core.sim1b[,c('r','m','mu_k','phi')]
true.param.1b <- true.param.1b[c('r','m','mu_k','phi')]

for (i in 1:4)
{
	for (j in 1:4)
	{
		if (j>i){
			plot(post.sample.core.sim1b[,i],post.sample.core.sim1b[,j],pch=1,col=rgb(0.67,0.84,0.9,0.1))
			HPDregionplot(mcmc(post.sample.core.sim1b[,c(i,j)]),add=TRUE,prob=c(0.5,0.95),n=100,lty=c(1,2))
			points(true.param.1b[i],true.param.1b[j],pch="+",col='red',cex=1.5)
		}
	}
}
dev.off()



# Figure S3 ----
load(here('results','post_jp_abot.RData'))
source(here('src','utility.R'))
pdf(file=here('figures_and_tables','figureS3.pdf'),width=5,height=4)
plot.fitted(r=post.sample.core.jp.abot[,'r'],m=post.sample.core.jp.abot[,'m'],mu_k=post.sample.core.jp.abot[,'mu_k'],timeRange=c(4000,1700),calendar = 'BCAD',xlab='',ylab='')
mtext('BCE/CE',side=1,line=2.2,cex=0.7)
mtext('Proportion',side=2,line=2.2,cex=0.7)
dev.off()

# Figure S4 ----
pdf(file=here('figures_and_tables','figureS4.pdf'),width=7,height=7)
layout(matrix(c(1,2,3,4,5,10,14,15,16,6,20,11,17,18,7,20,20,12,19,8,20,20,20,13,9),nrow=5,ncol=5),width=c(0.1,1,1,1,1),height=c(1,1,1,1,0.1))

label.plot <- function(label)
{
	par(mar=c(0,0,0,0))
	plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
	text(label,x=0.5,y=0.5,las=2,adj=0.5,cex=1.5)
}

label.plot('r')
label.plot('m')
label.plot(TeX(r'($\mu$)'))
label.plot(TeX(r'($\varphi$)'))

label.plot('')
label.plot('r')
label.plot('m')
label.plot(TeX(r'($\mu$)'))
label.plot(TeX(r'($\varphi$)'))

par(mar=c(2,2,2,2))
postHPDplot(post.sample.core.jp.abot$r,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
postHPDplot(post.sample.core.jp.abot$m,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
postHPDplot(post.sample.core.jp.abot$mu_k,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
postHPDplot(post.sample.core.jp.abot$phi,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)


post.sample.core.jp.abot  <- post.sample.core.jp.abot[,c('r','m','mu_k','phi')]


for (i in 1:4)
{
	for (j in 1:4)
	{
		if (j>i){
			plot(post.sample.core.jp.abot[,i],post.sample.core.jp.abot[,j],pch=1,col=rgb(0.67,0.84,0.9,0.1))
			HPDregionplot(mcmc(post.sample.core.jp.abot[,c(i,j)]),add=TRUE,prob=c(0.5,0.95),n=100,lty=c(1,2))
		}
	}
}
dev.off()




# Figure S5 ----
load(here('results','post_gb_abot.RData'))
source(here('src','utility.R'))
pdf(file=here('figures_and_tables','figureS5.pdf'),width=5,height=4)
plot.fitted(r=post.sample.core.gb.abot[,'r'],m=post.sample.core.gb.abot[,'m'],mu_k=post.sample.core.gb.abot[,'mu_k'],timeRange=c(7000,3000),calendar = 'BCAD',xlab='',ylab='')
mtext('BCE',side=1,line=2.2,cex=0.7)
mtext('Proportion',side=2,line=2.2,cex=0.7)
dev.off()

# Figure S6 ----

pdf(file=here('figures_and_tables','figureS6.pdf'),width=7,height=7)
layout(matrix(c(1,2,3,4,5,10,14,15,16,6,20,11,17,18,7,20,20,12,19,8,20,20,20,13,9),nrow=5,ncol=5),width=c(0.1,1,1,1,1),height=c(1,1,1,1,0.1))

label.plot <- function(label)
{
	par(mar=c(0,0,0,0))
	plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
	text(label,x=0.5,y=0.5,las=2,adj=0.5,cex=1.5)
}

label.plot('r')
label.plot('m')
label.plot(TeX(r'($\mu$)'))
label.plot(TeX(r'($\varphi$)'))

label.plot('')
label.plot('r')
label.plot('m')
label.plot(TeX(r'($\mu$)'))
label.plot(TeX(r'($\varphi$)'))

par(mar=c(2,2,2,2))
postHPDplot(post.sample.core.gb.abot$r,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
postHPDplot(post.sample.core.gb.abot$m,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
postHPDplot(post.sample.core.gb.abot$mu_k,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)
postHPDplot(post.sample.core.gb.abot$phi,xlab='',ylab='',axes=F)
axis(1,line=0,padj=-1,tck=-0.05,cex.axis=0.8)


post.sample.core.gb.abot  <- post.sample.core.gb.abot[,c('r','m','mu_k','phi')]


for (i in 1:4)
{
	for (j in 1:4)
	{
		if (j>i){
			plot(post.sample.core.gb.abot[,i],post.sample.core.gb.abot[,j],pch=1,col=rgb(0.67,0.84,0.9,0.1))
			HPDregionplot(mcmc(post.sample.core.gb.abot[,c(i,j)]),add=TRUE,prob=c(0.5,0.95),n=100,lty=c(1,2))
		}
	}
}
dev.off()


