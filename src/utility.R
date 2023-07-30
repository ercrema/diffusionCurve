logistic  <- function(x)
{
	p <- 1/(1+exp(-x))
	return(p)
}

sigmoid <- function(x,k,r,m)
{
	return(k/(1+exp(r*(x-m))))
} 

plot.fitted  <- function(r,m,mu_k,timeRange,calendar='BP',nsample=NULL,interval=0.95)
{
	# Setup
	require(rcarbon)
	if (length(unique(c(length(r),length(m),length(mu_k))))>1)
	{
		stop('Length of r, m, and mu_k should be the same')
	}
	if (is.null(nsample)){nsample  <- length(r)}
	i  <- sample(length(r),size=nsample)
	r.plot <- r[i]
	m.plot <- round(m[i])
	k.plot <- mu_k[i]

	# Fill Matrix
	mat  <- sapply(1:nsample,function(x,y,r,m,k){sigmoid(x=y,r=r[x],m=m[x],k=k[x])},y=timeRange[1]:timeRange[2],r=r.plot,m=m.plot,k=k.plot)

	# Calendar
	if (calendar=='BP')
	{
		ticks  <- rev(pretty(timeRange))
		ticks.loc  <- ticks
		calCaption  <- 'BP'
	}
	if (calendar=='BCAD')
	{
		ticks  <- pretty(BPtoBCAD(timeRange))
		if (any(ticks==0)){ticks[which(ticks==0)] <- 1}
		ticks.loc  <- BCADtoBP(ticks)
		ticks  <- abs(ticks)
		calCaption  <- 'BCE/CE'
	}
	mean.mat  <- apply(mat,1,mean)
	lo.mat  <- apply(mat,1,quantile,prob=(1-interval)/2)
	hi.mat  <- apply(mat,1,quantile,prob=interval + (1-interval)/2)
	plot(NULL,xlim=timeRange,ylim=c(0,1),axes=F,xlab=calCaption,ylab='Probability')
	polygon(c(timeRange[1]:timeRange[2],timeRange[2]:timeRange[1]),c(lo.mat,rev(hi.mat)),border=NA,col='lightblue')
	lines(timeRange[1]:timeRange[2],mean.mat,lwd=2,lty=2)
	axis(2)
	axis(1,at=ticks.loc,label=ticks)
	box()
}

ppcheck  <- function(x,obs,ppmat,timeRange)
{
	require(rcarbon)
	#x ... caldates
	#obs ... observed binary sequence
	#ppmat ... predicted binary sequence	
	spd.obs1  <- spd(x[which(obs==1)],timeRange=timeRange,verbose=F,spdnormalised = FALSE)
	spd.obs0  <- spd(x[which(obs==0)],timeRange=timeRange,verbose=F,spdnormalised = FALSE)
	plotyears  <- timeRange[1]:timeRange[2]
	obs.prop  <- spd.obs1[[2]][,2] / (spd.obs0[[2]][,2] + spd.obs1[[2]][,2]) 
	pred.prop  <- matrix(NA,ncol=ncol(ppmat),nrow=length(plotyears))
	pb <- txtProgressBar(min = 0, max = ncol(ppmat), style = 3, width = 50, char = "=")

	for (i in 1:ncol(ppmat))
	{
		setTxtProgressBar(pb, i)
		if (all(ppmat[,i]==0)){spd.sim1  <- rep(0,length(plotyears))} else {
		spd.sim1  <- spd(x[which(ppmat[,i]==1)],timeRange=timeRange,verbose=F,spdnormalised = F)[[2]][,2]}
		if (all(ppmat[,i]==1)){spd.sim0  <- rep(0,length(plotyears))} else {
		spd.sim0  <- spd(x[which(ppmat[,i]==0)],timeRange=timeRange,verbose=F,spdnormalised = F)[[2]][,2]}
		pred.prop[,i]  <- spd.sim1 / (spd.sim0 + spd.sim1) 
	}
	return(list(plotyears=plotyears,obs.prop=obs.prop,pred.prop=pred.prop))
}


plotPcheck <- function(x,calendar,interval=0.9,envelope.col='lightgrey',positive.col='red',negative.col='blue',obs.col='black',obs.lwd=2)
{
	require(rcarbon)
	if (calendar=='BP')
	{
		ticks  <- pretty(x$plotyears)
		tickLoc  <- ticks
		xlab  <- 'BP'
	}

	if (calendar=='BCAD')
	{

		ticks  <- pretty(BPtoBCAD(x$plotyears))
		if (any(ticks==0)){ticks[which(ticks==0)] <- 1}
		tickLoc  <- BCADtoBP(ticks)
		ticks  <- abs(ticks)
		xlab  <- 'BC/AD'
	}



	obs  <- x$obs.prop
	lo  <- apply(x$pred.prop,1,quantile,prob=(1-interval)/2,na.rm=T)
	hi  <- apply(x$pred.prop,1,quantile,prob= interval + (1-interval)/2,na.rm=T)
	NAcensors  <- which(is.na(obs)|is.na(lo)|is.na(hi))
	obs[NAcensors]  <- 0
	lo[NAcensors] <- 0
	hi[NAcensors]  <- 1


	# Boom and Bust Handling ####
	booms <- which(obs>hi)
	busts <- which(obs<lo)
	baseline <- rep(NA,length(obs))
	colpts = rep('grey',length(obs))
	colpts[booms] = 'red'
	colpts[busts] = 'blue'

	boomPlot <- baseline
	if (length(booms)>0){ boomPlot[booms]=obs[booms] }
	bustPlot <- baseline
	if (length(busts)>0){ bustPlot[busts]=obs[busts] }           

	boomBlocks <- vector("list")
	counter <- 0
	state <- "off"
	for (i in 1:length(boomPlot)){
		if (!is.na(boomPlot[i])&state=="off"){
			counter <- counter+1
			boomBlocks <- c(boomBlocks,vector("list",1))
			boomBlocks[[counter]] <- vector("list",2)
			boomBlocks[[counter]][[1]] <- boomPlot[i]
			boomBlocks[[counter]][[2]] <- x$plotyears[i]
			state <- "on"
		}
		if (state=="on"){
			if (!is.na(boomPlot[i])){
				boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[i])
				boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],x$plotyears[i])
			}
			if (is.na(boomPlot[i])){
				state <- "off"
			}
		}   
	}

	bustBlocks <- vector("list")
	counter <- 0
	state <- "off"
	for (i in 1:length(bustPlot)){
		if (!is.na(bustPlot[i])&state=="off"){
			counter <- counter+1
			bustBlocks <- c(bustBlocks,vector("list",1))
			bustBlocks[[counter]] <- vector("list",2)
			bustBlocks[[counter]][[1]] <- bustPlot[i]
			bustBlocks[[counter]][[2]] <- x$plotyears[i]
			state <- "on"
		}
		if (state=="on"){
			if (!is.na(bustPlot[i])){
				bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[i])
				bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],x$plotyears[i])
			}
			if (is.na(bustPlot[i])){
				state <- "off"
			}
		}   
	}

	plot(NULL,xlim=rev(range(x$plotyears)), ylim=c(0,1), type="n", col="white", ylab='Probability', xlab=xlab, xaxt="n")
	axis(1,at=tickLoc,labels=ticks)

	polygon(c(x$plotyears,rev(x$plotyears)),c(lo,rev(hi)),col=envelope.col,border=NA)

	if (length(booms)>0){
		for (i in 1:length(boomBlocks)){
			bbb = unique(boomBlocks[[i]][[2]])
			index = which(x$plotyears%in%bbb)
			polygon(c(bbb,rev(bbb)),c(obs[index],rev(hi[index])),border=NA,col=positive.col)
		}  
	}

	if (length(busts)>0){
		for (i in 1:length(bustBlocks)){
			bbb = unique(bustBlocks[[i]][[2]])
			index = which(x$plotyears%in%bbb)
			polygon(c(bbb,rev(bbb)),c(obs[index],rev(lo[index])),border=NA,col=negative.col)
		}  
	}
	lines(x$plotyears,obs,lwd=obs.lwd,col=obs.col)
	if (length(NAcensors)>0)
	{
		Breaks <- c(0, which(diff(NAcensors) != 1), length(NAcensors))
		NAblocks <- sapply(seq(length(Breaks) - 1), function(i) NAcensors[c((Breaks[i] + 1),Breaks[i+1])])
		for(i in 1:ncol(NAblocks))
		{
			rect(xleft=x$plotyears[NAblocks[1,i]],xright=x$plotyears[NAblocks[2,i]],ybottom=-2,ytop=2,col='darkgrey',border='darkgrey')
		}
	}
	box()
}
