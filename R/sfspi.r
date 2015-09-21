#pi from an SFS assumed to be log-scaled

sfspi<-function(x,fold=TRUE){
	if(fold){
		nchr<-(length(x)-1)*2
		nbin<-length(x)
		nind<-length(x)-1
		denom<-(nchr*(nchr-1))/2
		pi<-sum(((1:nind)*(nchr-(1:nind)))*exp(x[2:nbin]))/denom
		return(pi)
		}
	if(!fold){
		nchr<-length(x)-1
		nbin<-length(x)
		nind<-nchr/2
		denom<-(nchr*(nchr-1))/2
		pi<-sum(((0:nchr)*(nchr-(0:nchr)))*exp(x))/denom
		return(pi)
		}
	}
