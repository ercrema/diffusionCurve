logistic  <- function(x)
{
	p <- 1/(1+exp(-x))
	return(p)
}

sigmoid <- function(x,k,r,m)
{
	return(k/(1+exp(r*(x-m))))
} 

