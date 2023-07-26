# Figure 1 ----
library(here)
library(rnaturalearth)
library(sf)
library(dplyr)
library(RColorBrewer)
cols  <- brewer.pal('Set2',n=6)

load(here('data','gbdata.RData'))
load(here('data','burialdata.RData'))
load(here('data','ricedata.RData'))

win  <- ne_countries(scale=10,returnclass='sf') |> st_combine()

japan  <- ne_states(country = 'japan',returnclass = 'sf') |> subset(!name%in%c('Hokkaidō','Okinawa'))
japan$Region  <- 'SW Japan'
japan$Region[which(japan$name%in%c('Aomori','Akita','Iwate','Miyagi','Yamagata','Fukushima'))]  <- 'NE Japan'
japan$Region[which(japan$name%in%c('Fukui','Ishikawa','Toyama','Gifu','Aichi','Shizuoka','Nagano','Niigata','Yamanashi','Kanagawa','Tokyo','Ibaraki','Gunma','Tochigi','Chiba','Saitama'))]  <- 'C Japan'
japan  <- group_by(japan,by=Region) |> summarise()
japan.sites  <- select(ricedata,Latitude,Longitude) |> na.exclude() |> unique() |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)

uk <- ne_states(country = "united kingdom",returnclass = 'sf')
scotland  <- subset(uk,region%in%c('Highlands and Islands','North Eastern','South Western',"Eastern")) |> st_union() |> st_cast(to='POLYGON')
scotland_mainland  <- scotland[order(st_area(scotland),decreasing=TRUE)[1]] |> st_sf()
scottish_isles  <- scotland[order(st_area(scotland),decreasing=TRUE)[-1]] |> st_sf()
england_and_wales  <- subset(uk,!region%in%c('Highlands and Islands','North Eastern','South Western', 'Eastern','Northern Ireland')) |> st_union() |> st_cast(to='POLYGON') |> st_sf()
gb  <- st_union(scotland_mainland,scottish_isles,england_and_wales)
gb.abot.sites  <- select(gbdata,Longitude,Latitude) |> na.exclude() |> unique() |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)
gb.burial.sites  <- select(burial,Longitude,Latitude) |> na.exclude() |> unique() |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)



pdf(file = here('figures','figure1.pdf'),width=5.5,height=3.5,pointsize=1)
par(mfrow=c(1,2)) 
plot(win,xlim=c(128,143),ylim=c(33,41),lwd=0.5,border='darkgrey',col='darkgrey',bg='powderblue')
abline(v=seq(120,144,2),h=seq(26,44,2),col='grey77',lwd=0.5)
plot(japan,col=cols[1:3],add=T,border=NA)
plot(japan.sites,pch=20,add=T,cex=0.8)
axis(1,tck=-0.03,padj=-0.5,at=seq(120,144,4))
axis(1,tck=-0.01,at=seq(120,144,2),labels=NA)
axis(2,tck=-0.03,padj=-0.5,at=seq(26,44,4))
axis(2,tck=-0.02,at=seq(26,44,2),labels=NA)
mtext('Longitide',side=1,line=2.5)
mtext('Latitude',side=2,line=2.5)
box()
legend(x=132,y=41.5,legend=c('SW Japan','C Japan','NE Japan'),fill=cols[c(3,1,2)],bty='n',cex=.9)
legend(x=132,y=39.5,legend=c('Sampling Sites'),pch=20,bty='n',cex=.9)

plot(win,xlim=c(-6,1),ylim=c(50,61),border='darkgrey',col='darkgrey',bg='powderblue')
abline(v=seq(-10,4,2),h=seq(48,64,2),col='grey77',lwd=0.5)
plot(scotland_mainland,add=TRUE,col=cols[5],border=NA)
plot(scottish_isles,add=TRUE,col=cols[6],border=NA)
plot(england_and_wales,add=TRUE,col=cols[4],border=NA)
plot(gb.abot.sites,pch=20,add=T,cex=0.8)
plot(gb.burial.sites,pch=3,add=T,cex=0.8)
axis(2,tck=-0.03,padj=-0.5,at=seq(48,64,4))
axis(2,tck=-0.01,at=seq(48,64,2),labels=NA)
axis(1,tck=-0.03,padj=-0.5,at=seq(-12,12,4))
axis(1,tck=-0.02,at=seq(-12,12,2),labels=NA)
mtext('Longitide',side=1,line=2.5)
mtext('Latitude',side=2,line=2.5)
box()
legend(x=-1.5,y=58,legend=c('Scottish Isles','Scotland Mainland','England & Wales'),fill=cols[6:4],bty='n',cex=.9)
legend(x=-1.5,y=56.5,legend=c('Sampling Sites (Seeds)','Sampling Sites (Burials)'),pch=c(20,3),bty='n',cex=.9)
dev.off()


# Figure 2 (Fitted Parameters GB) ----
library(here)
library(rcarbon)
library(ggplot2)
library(tidybayes)
library(tidyr)
library(gridExtra)
load(here('results','post_gbwhole.RData'))
load(here('results','post_gbregional.RData'))
source(here('src','utility.R'))

# Posterior r
post.r.gb  <- post.sample.core.gbregional[,grep('r\\[',colnames(post.sample.core.gbregional))]
post.r.gb  <- cbind(post.r.gb,post.sample.core.gbwhole[,'r'])
colnames(post.r.gb)  <- c('England and Wales','Scotland (mainland)','Scottish Isles','Britain (whole)')
post.r.gb  <- gather(post.r.gb)
r.britain <- ggplot(aes(x = key, y = value),data=post.r.gb) + 
	     stat_interval(aes(x = key), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	     scale_color_brewer() +
	     xlab('')   +
	     ylab('Adoption Rate (r)') +
	     ggtitle('a') +
	     theme(axis.text.y = element_text(angle = 90))

# Posterior mu_k
post.k.gb  <- post.sample.core.gbregional[,grep('mu_k\\[',colnames(post.sample.core.gbregional))]
post.k.gb  <- cbind(post.k.gb,post.sample.core.gbwhole[,'mu_k'])
colnames(post.k.gb)  <- c('England and Wales','Scotland (mainland)','Scottish Isles','Britain (whole)')
post.k.gb  <- gather(post.k.gb)
post.k.gb$value  <- logistic(post.k.gb$value)

k.britain  <- ggplot(aes(x = key, y = value),data=post.k.gb) + 
	      stat_interval(aes(x = key), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	      scale_color_brewer() +
	      xlab('')   +
	      ylab('Adoption Cap (average k)') +
	      ggtitle('b') +
	      theme(axis.text.y = element_text(angle = 90))

# Posterior m
post.m.gb  <- post.sample.core.gbregional[,grep('m\\[',colnames(post.sample.core.gbregional))]
post.m.gb  <- cbind(post.m.gb,post.sample.core.gbwhole[,'m'])
colnames(post.m.gb)  <- c('England and Wales','Scotland (mainland)','Scottish Isles','Britain (whole)')
post.m.gb  <- gather(post.m.gb) 
post.m.gb$value = BPtoBCAD(round(post.m.gb$value)) #negative as tidybayes struggles in reversing x-axis

# pdf(here('figures','post_m_britain.pdf'),width=7,height=6)
m.britain <- ggplot(aes(y = value, x = key),data=post.m.gb) + 
	     stat_interval(aes(y = value), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	     scale_color_brewer() +
	     xlab('')   +
	     ylab('Midpoint (m) in BP') +
	     ggtitle('c') +
	     scale_y_continuous(labels = abs)

grid.arrange(r.britain,k.britain,m.britain,ncol=1)


# Figure 2 (Fitted Parameters Japan) ----
library(here)
library(rcarbon)
library(ggplot2)
library(tidybayes)
library(tidyr)
library(gridExtra)
load(here('results','post_jpwhole.RData'))
load(here('results','post_jpregional.RData'))
source(here('src','utility.R'))

# Posterior r
post.r.jp  <- post.sample.core.jpregional[,grep('r\\[',colnames(post.sample.core.jpregional))]
post.r.jp  <- cbind(post.r.jp,post.sample.core.jpwhole[,'r'])
colnames(post.r.jp)  <- c('SW Japan','C Japan','NE Japan','Japan (excluding Hokkaido)')
post.r.jp  <- gather(post.r.jp)
r.japan <- ggplot(aes(x =key , y = value),data=post.r.jp) + 
	     stat_interval(aes(x =factor(key,level=c('Japan (excluding Hokkaido)','SW Japan','C Japan','NE Japan')), ), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	     scale_color_brewer() +
	     xlab('')   +
	     ylab('Adoption Rate (r)') +
	     ggtitle('a') +
	     theme(axis.text.y = element_text(angle = 90))

# Posterior mu_k
post.k.jp  <- post.sample.core.jpregional[,grep('mu_k\\[',colnames(post.sample.core.jpregional))]
post.k.jp  <- cbind(post.k.jp,post.sample.core.jpwhole[,'mu_k'])
colnames(post.k.jp)  <- c('SW Japan','C Japan','NE Japan','Japan (excluding Hokkaido)')
post.k.jp  <- gather(post.k.jp)
post.k.jp$value  <- logistic(post.k.jp$value)

k.japan  <- ggplot(aes(x = key, y = value),data=post.k.jp) + 
	      stat_interval(aes(x = factor(key,levels=c('Japan (excluding Hokkaido)','SW Japan','C Japan','NE Japan'))), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	      scale_color_brewer() +
	      xlab('')   +
	      ylab('Adoption Cap (average k)') +
	      ggtitle('b') +
	      theme(axis.text.y = element_text(angle = 90))

# Posterior m
post.m.jp  <- post.sample.core.jpregional[,grep('m\\[',colnames(post.sample.core.jpregional))]
post.m.jp  <- cbind(post.m.jp,post.sample.core.jpwhole[,'m'])
colnames(post.m.jp)  <- c('SW Japan','C Japan','NE Japan','Japan (excluding Hokkaido)')
post.m.jp  <- gather(post.m.jp) 
post.m.jp$value = BPtoBCAD(round(post.m.jp$value)) #negative as tidybayes struggles in reversing x-axis

# pdf(here('figures','post_m_japan.pdf'),width=7,height=6)
m.japan <- ggplot(aes(y = value, x = factor(key,levels=c('Japan (excluding Hokkaido)','SW Japan','C Japan','NE Japan'))),data=post.m.jp) + 
	     stat_interval(aes(y = value), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	     scale_color_brewer() +
	     xlab('')   +
	     ylab('Midpoint (m) in BP') +
	     ggtitle('c') +
	     scale_y_continuous(labels = abs)

grid.arrange(r.japan,k.japan,m.japan,ncol=1)

# Figure 2 (Fitted Parameters Japan v2) ----
library(here)
library(rcarbon)
library(ggplot2)
library(tidybayes)
library(tidyr)
library(gridExtra)
load(here('results','post_jpwhole.RData'))
load(here('results','post_jpregional.RData'))
source(here('src','utility.R'))

# Posterior r
post.r.jp  <- post.sample.core.jpregional[,grep('r\\[',colnames(post.sample.core.jpregional))]
post.r.jp  <- cbind(post.r.jp,post.sample.core.jpwhole[,'r'])
colnames(post.r.jp)  <- c('SW Japan','C Japan','NE Japan','Japan (excluding Hokkaido)')
post.r.jp  <- gather(post.r.jp)
r.japan <- ggplot(aes(x =key , y = value),data=post.r.jp) + 
	     stat_interval(aes(x =factor(key,level=c('Japan (excluding Hokkaido)','SW Japan','C Japan','NE Japan')), ), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	     scale_color_brewer() +
	     xlab('')   +
	     ylab('Adoption Rate (r)') +
	     ggtitle('a') +
	     theme(axis.text.y = element_text(angle = 90))

# Posterior mu_k
post.k.jp  <- post.sample.core.jpregional[,grep('mu_k\\[',colnames(post.sample.core.jpregional))]
post.k.jp  <- cbind(post.k.jp,post.sample.core.jpwhole[,'mu_k'])
colnames(post.k.jp)  <- c('SW Japan','C Japan','NE Japan','Japan (excluding Hokkaido)')
post.k.jp  <- gather(post.k.jp)
post.k.jp$value  <- logistic(post.k.jp$value)

k.japan  <- ggplot(aes(x = key, y = value),data=post.k.jp) + 
	      stat_interval(aes(x = factor(key,levels=c('Japan (excluding Hokkaido)','SW Japan','C Japan','NE Japan'))), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	      scale_color_brewer() +
	      xlab('')   +
	      ylab('Adoption Cap (average k)') +
	      ggtitle('b') +
	      theme(axis.text.y = element_text(angle = 90))

# Posterior m
post.m.jp  <- post.sample.core.jpregional[,grep('m\\[',colnames(post.sample.core.jpregional))]
post.m.jp  <- cbind(post.m.jp,post.sample.core.jpwhole[,'m'])
colnames(post.m.jp)  <- c('SW Japan','C Japan','NE Japan','Japan (excluding Hokkaido)')
post.m.jp  <- gather(post.m.jp) 
post.m.jp$value = BPtoBCAD(round(post.m.jp$value)) #negative as tidybayes struggles in reversing x-axis

# pdf(here('figures','post_m_japan.pdf'),width=7,height=6)
m.japan <- ggplot(aes(y = value, x = factor(key,levels=c('Japan (excluding Hokkaido)','SW Japan','C Japan','NE Japan'))),data=post.m.jp) + 
	     stat_interval(aes(y = value), .width = c(.5, .8, .95), interval_size_range = c(5,8)) +
	     scale_color_brewer() +
	     xlab('')   +
	     ylab('Midpoint (m) in BP') +
	     ggtitle('c') +
	     scale_y_continuous(labels = abs)

grid.arrange(r.japan,k.japan,m.japan,ncol=1)
# Figure 3 ----

load(here('results','ppc_jpwhole.RData'))
load(here('data','ricedata.RData'))
cal <- calibrate(ricedata$C14Age,ricedata$C14Error) |> medCal()

plotPcheck(ppc.jpwhole,calendar='BCAD')
points(x=cal[which(ricedata$rice_nuts==1)],y=rep(0.98,times=sum(ricedata$rice_nuts==1)),pch="+",col=rgb(0,0,0,0.2))
points(x=cal[which(ricedata$rice_nuts==0)],y=rep(0.02,times=sum(ricedata$rice_nuts==0)),pch="+",col=rgb(0,0,0,0.2))



# Figure 4 (ICAR Model)

