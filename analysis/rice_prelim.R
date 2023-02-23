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


ricedata$Region3  <- 'SW'
ricedata$Region3[which(ricedata$Region %in% c('Chubu','Kanto'))] = 'C'
ricedata$Region3[which(ricedata$Region %in% c('Tohoku'))] = 'NE'

ricedata.sw3  <- subset(ricedata,Region3=='SW')
ricedata.c3  <- subset(ricedata,Region3=='C')
ricedata.ne3  <- subset(ricedata,Region3=='NE')

caldates.sw3 <- calibrate(ricedata.sw3$C14Age,ricedata.sw3$C14Error) 
caldates.c3 <- calibrate(ricedata.c3$C14Age,ricedata.c3$C14Error)
caldates.ne3 <- calibrate(ricedata.ne3$C14Age,ricedata.ne3$C14Error)

stspd.sw3  <- stackspd(caldates.sw3,runm = 20,timeRange=c(4000,1750),group = ricedata.sw3$group)
stspd.ne3  <- stackspd(caldates.ne3,runm = 20,timeRange=c(4000,1750),group = ricedata.ne3$group)
stspd.c3  <- stackspd(caldates.c3,runm = 20,timeRange=c(4000,1750),group = ricedata.c3$group)


par(mfrow=c(3,1))
plot(stspd.sw3,type='proportion',col.fill=c('steelblue','lightgrey'),main='SW Japan')
plot(stspd.c3,type='proportion',col.fill=c('steelblue','lightgrey'),main='C Japan')
plot(stspd.ne3,type='proportion',col.fill=c('steelblue','lightgrey'),main='NE Japan')
