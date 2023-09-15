# Load Libraries ----
library(here)
library(rnaturalearth)
library(sf)
library(dplyr)
library(rcarbon)
library(nimbleCarbon)
library(RColorBrewer)

# Load R images and custom functions ----
source(here('src','utility.R'))
load(here('data','gbdata.RData'))
load(here('data','burialdata.RData'))
load(here('data','jpdata.RData'))
load(here('results','ppc_jp_abot.RData'))
load(here('results','post_jp_abot.RData'))
load(here('results','ppc_gb_abot.RData'))
load(here('results','post_gb_abot.RData'))
load(here('results','post_icar_burial.RData'))
load(here('sim','simdata','simdata1.RData'))
load(here('sim','simdata','simdata2.RData'))
load(here('sim','simdata','simdata3.RData'))
load(here('sim','results','post_sim1.RData'))
load(here('sim','results','post_sim2.RData'))
load(here('sim','results','post_icar_sim3.RData'))


# Figure 1 (Distribution Maps) ----
win  <- ne_countries(scale=10,returnclass='sf') |> st_combine()
japan  <- ne_states(country = 'japan',returnclass = 'sf') |> subset(!name%in%c('Hokkaidō','Okinawa')) |> st_union()
japan.sites  <- select(jpdata,Latitude,Longitude) |> na.exclude() |> unique() |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)

gb <- ne_states(country = "united kingdom",returnclass = 'sf') |> subset(region!='Northern Ireland') |> st_union()
gb.abot.sites  <- select(gbdata,Longitude,Latitude) |> na.exclude() |> unique() |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)
gb.burial.sites  <- select(burial,Longitude,Latitude) |> na.exclude() |> unique() |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)

pdf(file = here('figures_and_tables','figure1.pdf'),width=5.5,height=3.5,pointsize=1)
par(mfrow=c(1,2)) 
plot(win,xlim=c(128,143),ylim=c(33,41),lwd=0.5,border='darkgrey',col='darkgrey',bg='powderblue')
abline(v=seq(120,144,2),h=seq(26,44,2),col='grey77',lwd=0.5)
plot(japan,col='cornsilk',add=T,border=NA)
plot(japan.sites,pch=1,add=T,cex=0.8)
axis(1,tck=-0.03,padj=-0.5,at=seq(120,144,4))
axis(1,tck=-0.01,at=seq(120,144,2),labels=NA)
axis(2,tck=-0.03,padj=-0.5,at=seq(26,44,4))
axis(2,tck=-0.02,at=seq(26,44,2),labels=NA)
mtext('Longitide',side=1,line=2.5)
mtext('Latitude',side=2,line=2.5)
box()
legend('topleft',legend='Sampling Sites',pch=1,bg='white')

plot(win,xlim=c(-8,1),ylim=c(50,61),border='darkgrey',col='darkgrey',bg='powderblue')
abline(v=seq(-10,4,2),h=seq(48,64,2),col='grey77',lwd=0.5)
plot(gb,col='cornsilk',add=T,border=NA)
plot(gb.abot.sites,pch=1,add=T,cex=0.8)
plot(gb.burial.sites,pch=2,add=T,cex=0.8)
axis(2,tck=-0.03,padj=-0.5,at=seq(48,64,4))
axis(2,tck=-0.01,at=seq(48,64,2),labels=NA)
axis(1,tck=-0.03,padj=-0.5,at=seq(-12,12,4))
axis(1,tck=-0.02,at=seq(-12,12,2),labels=NA)
mtext('Longitide',side=1,line=2.5)
mtext('Latitude',side=2,line=2.5)
box()
legend(x='topleft',legend=c('Sampling Sites (Seeds)','Sampling Sites (Burials)'),pch=c(1,2),bg='white')
dev.off()

# Figure 2 (Fitted model, simulations 1-3) ----
pdf(file = here('figures_and_tables','figure2.pdf'),width=3.5,height=4.5,pointsize=9)
par(mfrow=c(3,1),mar=c(4,3.5,0.3,0.1),lend=2)

plot.fitted(r=post.sample.core.sim1[,'r'],m=post.sample.core.sim1[,'m'],mu_k=post.sample.core.sim1[,'mu_k'],timeRange=c(4000,1700),calendar = 'BCAD',xlab='',ylab='')
mtext('BCE/CE',side=1,line=2.2,cex=0.7)
mtext('Proportion',side=2,line=2.2,cex=0.7)
lines(4000:1700,sigmoid(x=4000:1700,m=true.param.1$m,k=true.param.1$mu_k,r=true.param.1$r),lty=3,col=2,lwd=2)
legend('topright',legend=c('a'),bty='n',cex=1.2)


plot.fitted(r=post.sample.core.sim2[,'r'],m=post.sample.core.sim2[,'m'],mu_k=post.sample.core.sim2[,'mu_k'],timeRange=c(7000,3000),calendar = 'BCAD',xlab='',ylab='')
mtext('BCE',side=1,line=2.2,cex=0.7)
mtext('Proportion',side=2,line=2.2,cex=0.7)
lines(7000:3000,sigmoid(x=7000:3000,m=true.param.2$m,k=true.param.2$mu_k,r=true.param.2$r),lty=3,col=2,lwd=2)
legend(x=7000,y=1,legend=c('Model','Posterior Mean','95% HPD'),lty=c(3,2,1),col=c(2,1,'lightblue'),lwd=c(1.5,1.5,4),bty='n')
legend('topright',legend=c('b'),bty='n',cex=1.2)


tblocks.sim3 <- seq(constants.icar.sim3$a,constants.icar.sim3$b,by=-constants.icar.sim3$res)
plot(tblocks.sim3,apply(post.sample.combined.icar.sim3[,1:constants.icar.sim3$n.tblocks],2,mean),pch=20,xlim=c(constants.icar.sim3$a,constants.icar.sim3$b),ylim=c(0,1),type='n',xlab='',ylab='',axes=F)
mtext('BCE',side=1,line=2.2,cex=0.7)
mtext('Proportion',side=2,line=2.2,cex=0.7)
for (i in 1:constants.icar.sim3$n.tblocks)
{
	rect(xleft=tblocks.sim3[i]+45,xright=tblocks.sim3[i]-45,ybottom=quantile(post.sample.combined.icar.sim3[,i],0.025),ytop=quantile(post.sample.combined.icar.sim3[,i],0.975),border=NA,col='lightblue')
	rect(xleft=tblocks.sim3[i]+45,xright=tblocks.sim3[i]-45,ybottom=quantile(post.sample.combined.icar.sim3[,i],0.25),ytop=quantile(post.sample.combined.icar.sim3[,i],0.75),border=NA,col='steelblue')

}
lines(tblocks.sim3,apply(post.sample.combined.icar.sim3[,1:constants.icar.sim3$n.tblocks],2,mean),pch=20,type='b')
lines(timeRange[1]:timeRange[2],Pseq,col=2,lty=3,lwd=2)
legend(x=4950,y=1.05,legend=c('Model','Posterior Mean','95% HPD','50% HPD'),pch=c(NA,20,NA,NA),lty=c(3,1,1,1),col=c(2,1,'lightblue','steelblue'),lwd=c(1.5,1.5,4,4),bty='n')
legend('topright',legend=c('c'),bty='n',cex=1.2)
box()
axis(1)
axis(2)
dev.off()



# Figure 3 (Posterior Predictive Checks/Japan) ----
cal.jp <- calibrate(jpdata$C14Age,jpdata$C14Error) |> medCal()
pdf(file = here('figures_and_tables','figure3.pdf'),width=5.5,height=4,pointsize=8)
plotPcheck(ppc.jp.abot,calendar='BCAD')
axis(3)
axis(3,at=c(seq(5000,100,-100)),labels=NA,tck=-0.01)
axis(1,at=BCADtoBP(c(seq(-3000,-100,100),seq(100,1000,100))),labels=NA,tck=-0.01)
mtext(side=3,line=3,'cal BP')
points(x=cal.jp[which(jpdata$rice_nuts==1)],y=rep(0.98,times=sum(jpdata$rice_nuts==1)),pch="+",col=rgb(0,0,0,0.2))
points(x=cal.jp[which(jpdata$rice_nuts==0)],y=rep(0.02,times=sum(jpdata$rice_nuts==0)),pch="+",col=rgb(0,0,0,0.2))
legend(x=4000,y=0.95,legend=c('Observed SPD','Posterior Mean'),lty=c(1,2),lwd=c(2,1),bty='n')
legend(x=4000,y=0.85,legend=c('90% Prediction Envelope','Positive Deviation','Negative Deviation'),fill=c('lightgrey','red','blue'),bty='n')
legend(x=3972,y=0.68,legend=c('Median Dates'),pch=3,bty='n',cex=1)
dev.off()

# Figure 4 (Posterior Predictive Checks/Britain) ----
cal.gb <- calibrate(gbdata$CRA,gbdata$Error) |> medCal()
pdf(file = here('figures_and_tables','figure4.pdf'),width=5.5,height=4,pointsize=8)
plotPcheck(ppc.gb.abot,calendar='BCAD')
axis(3)
axis(3,at=c(seq(8000,100,-100)),labels=NA,tck=-0.01)
axis(1,at=BCADtoBP(c(seq(-6000,-100,100),seq(100,1000,100))),labels=NA,tck=-0.01)
mtext(side=3,line=3,'cal BP')
points(x=cal.gb[which(gbdata$cat2=='Wheat + Barley')],y=rep(0.98,times=sum(gbdata$cat2=='Wheat + Barley')),pch="+",col=rgb(0,0,0,0.2))
points(x=cal.gb[which(gbdata$cat2=='Hazelnut')],y=rep(0.02,times=sum(gbdata$cat2=='Hazelnut')),pch="+",col=rgb(0,0,0,0.2))
legend(x=7000,y=0.95,legend=c('Observed SPD','Posterior Mean'),lty=c(1,2),lwd=c(2,1),bty='n')
legend(x=7000,y=0.85,legend=c('90% Prediction Envelope','Positive Deviation','Negative Deviation'),fill=c('lightgrey','red','blue'),bty='n')
legend(x=6972,y=0.68,legend=c('Median Dates'),pch=3,bty='n',cex=1)
dev.off()

# Figure 5 (Posterior ICAR) ----
tblocks <- seq(constants.icar$a,constants.icar$b,by=-constants.icar$res)
cal.burial <- calibrate(burial$CRA,burial$Error) |> medCal()


pdf(here('figures_and_tables','figure5.pdf'),width=5.5,height=4.5,pointsize=8)
plot(tblocks,apply(post.sample.combined.icar[,1:constants.icar$n.tblocks],2,mean),pch=20,xlim=c(constants.icar$a,constants.icar$b),ylim=c(0,1),type='n',xlab='BCE',ylab='Estimated Proportion of Cremations',axes=F)
for (i in 1:constants.icar$n.tblocks)
{
	rect(xleft=tblocks[i]+45,xright=tblocks[i]-45,ybottom=quantile(post.sample.combined.icar[,i],0.025),ytop=quantile(post.sample.combined.icar[,i],0.975),border=NA,col='lightblue')
	rect(xleft=tblocks[i]+45,xright=tblocks[i]-45,ybottom=quantile(post.sample.combined.icar[,i],0.25),ytop=quantile(post.sample.combined.icar[,i],0.75),border=NA,col='steelblue')

}
lines(tblocks,apply(post.sample.combined.icar[,1:constants.icar$n.tblocks],2,mean),pch=20,type='b')
legend(x=2600,y=0.85,legend=c('95% HPD','50% HPD'),fill=c('lightblue','steelblue',1),cex=1,bty='n')
legend(x=2600,y=0.75,legend=c('Median Dates'),pch='+',bty='n',cex=1)
axis(3)
mtext(side=3,line=3,'cal BP')
axis(1,at=BCADtoBP(c(seq(-3500,-500,500),1)),labels=c(seq(3500,500,-500),1))
axis(2)
axis(3,at=c(seq(8000,100,-100)),labels=NA,tck=-0.01)
axis(1,at=BCADtoBP(c(seq(-6000,-100,100),seq(100,1000,100))),labels=NA,tck=-0.01)
box()
points(x=cal.burial[which(burial$Cremation==1)],y=rep(1,times=sum(burial$Cremation==1)),pch="+",col=rgb(0,0,0,0.2))
points(x=cal.burial[which(burial$Cremation==0)],y=rep(0,times=sum(burial$Cremation==0)),pch="+",col=rgb(0,0,0,0.2))
dev.off()

