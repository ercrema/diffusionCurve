# Figure 1 ----
library(here)
library(rnaturalearth)
library(sf)
library(dplyr)
library(RColorBrewer)

load(here('data','gbdata.RData'))
load(here('data','burialdata.RData'))
load(here('data','jpdata.RData'))

win  <- ne_countries(scale=10,returnclass='sf') |> st_combine()

japan  <- ne_states(country = 'japan',returnclass = 'sf') |> subset(!name%in%c('Hokkaidō','Okinawa')) |> st_union()
japan.sites  <- select(jpdata,Latitude,Longitude) |> na.exclude() |> unique() |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)

uk <- ne_states(country = "united kingdom",returnclass = 'sf') |> subset(region!='Northern Ireland') |> st_union()
gb.abot.sites  <- select(gbdata,Longitude,Latitude) |> na.exclude() |> unique() |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)
gb.burial.sites  <- select(burial,Longitude,Latitude) |> na.exclude() |> unique() |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)



pdf(file = here('figures','figure1.pdf'),width=5.5,height=3.5,pointsize=1)
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
plot(uk,col='cornsilk',add=T,border=NA)
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

# Figure 2 (Posterior Japan) ----
library(nimbleCarbon)
load(here('results','post_jp_abot.RData'))
r.jp  <- post.sample.core.jp.abot[,'r']
k.jp  <- post.sample.core.jp.abot[,'mu_k'] |> logistic()
m.jp  <- post.sample.core.jp.abot[,'m'] |> round()

# Figure 3 (Posterior Predictive Checks/Japan) ----
library(here)
library(rcarbon)
source(here('src','utility.R'))
load(here('results','ppc_jp_abot.RData'))
load(here('results','post_jp_abot.RData'))
load(here('data','jpdata.RData'))
cal.jp <- calibrate(jpdata$C14Age,jpdata$C14Error) |> medCal()
r.jp  <- post.sample.core.jp.abot[,'r']
k.jp  <- post.sample.core.jp.abot[,'mu_k'] |> logistic()
m.jp  <- post.sample.core.jp.abot[,'m'] |> round()
mat.jp  <- sapply(1:nrow(post.sample.core.jp.abot),function(x,y,r,m,k){sigmoid(x=y,r=r[x],m=m[x],k=k[x])},y=4000:1700,r=r.jp,m=m.jp,k=k.jp)
pdf(file = here('figures','figure3.pdf'),width=5.5,height=4,pointsize=8)
plotPcheck(ppc.jp.abot,calendar='BCAD')
lines(4000:1700,apply(mat.jp,1,mean),lty=2)
points(x=cal.jp[which(jpdata$rice_nuts==1)],y=rep(0.98,times=sum(jpdata$rice_nuts==1)),pch="+",col=rgb(0,0,0,0.2))
points(x=cal.jp[which(jpdata$rice_nuts==0)],y=rep(0.02,times=sum(jpdata$rice_nuts==0)),pch="+",col=rgb(0,0,0,0.2))
dev.off()



# Figure 4 (Posterior Britain) ----
library(nimbleCarbon)

# Figure 5 (Posterior Predictive Checks/Britain) ----
library(here)
library(rcarbon)
source(here('src','utility.R'))
load(here('results','ppc_gb_abot.RData'))
load(here('results','post_gb_abot.RData'))
load(here('data','gbdata.RData'))
cal.gb <- calibrate(gbdata$CRA,gbdata$Error) |> medCal()
r.gb  <- post.sample.core.gb.abot[,'r']
k.gb  <- post.sample.core.gb.abot[,'mu_k'] |> logistic()
m.gb  <- post.sample.core.gb.abot[,'m'] |> round()
mat.gb  <- sapply(1:nrow(post.sample.core.gb.abot),function(x,y,r,m,k){sigmoid(x=y,r=r[x],m=m[x],k=k[x])},y=7000:3000,r=r.gb,m=m.gb,k=k.gb)

pdf(file = here('figures','figure5.pdf'),width=5.5,height=4,pointsize=8)
plotPcheck(ppc.gb.abot,calendar='BCAD')
lines(7000:3000,apply(mat.gb,1,mean),lty=2)
points(x=cal.gb[which(gbdata$cat2=='Wheat + Barley')],y=rep(0.98,times=sum(gbdata$cat2=='Wheat + Barley')),pch="+",col=rgb(0,0,0,0.2))
points(x=cal.gb[which(gbdata$cat2=='Hazelnut')],y=rep(0.02,times=sum(gbdata$cat2=='Hazelnut')),pch="+",col=rgb(0,0,0,0.2))
dev.off()



