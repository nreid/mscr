#this function should fold a log-scaled sfs

foldsfs<-function(x,islog=TRUE){
	if(!islog){stop("sfs must be log-scaled")}
  nchr<-length(x)-1
	nbin<-length(x)
	nind<-nchr/2

	x2<-exp(x)
	x3<-x2[1:(nind+1)]+c(x2[nbin:(nind+2)],0)
	x3<-log(x3)
	return(x3)
	}

