# Load library
library(here)
library(rcarbon)
load(here('data','ricedata.RData'))
ricedata  <- ricedata[order(ricedata$rice_nuts,decreasing = TRUE),]
ricedata$group  <- ifelse(ricedata$rice_nuts==1,'rice','nuts')
ricedata$group  <- factor(ricedata$group,levels=c('rice','nuts'),ordered=TRUE)
caldates <- calibrate(ricedata$C14Age,ricedata$C14Error)
stspd  <- stackspd(caldates,runm = 20,timeRange=c(4000,1750),group = ricedata$group)
plot(stspd,type='proportion',col.fill=c('steelblue','lightgrey'))

ricedata$Region2  <- 'SW'
ricedata$Region2[which(ricedata$Region %in% c('Chubu','Kanto','Tohoku'))] = 'NE'
ricedata.sw  <- subset(ricedata,Region2=='SW')
ricedata.ne  <- subset(ricedata,Region2=='NE')

caldates.sw <- calibrate(ricedata.sw$C14Age,ricedata.sw$C14Error) 
caldates.ne <- calibrate(ricedata.ne$C14Age,ricedata.ne$C14Error)
stspd.sw  <- stackspd(caldates.sw,runm = 20,timeRange=c(4000,1750),group = ricedata.sw$group)
stspd.ne  <- stackspd(caldates.ne,runm = 20,timeRange=c(4000,1750),group = ricedata.ne$group)

par(mfrow=c(2,1))
plot(stspd.sw,type='proportion',col.fill=c('steelblue','lightgrey'),main='SW Japan')
plot(stspd.ne,type='proportion',col.fill=c('steelblue','lightgrey'),main='NE Japan')

