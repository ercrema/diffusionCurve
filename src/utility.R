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
	k.plot <- logistic(mu_k[i])

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








