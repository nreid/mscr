#watterson's theta from a log-scaled SFS

sfsw<-function(x,fold=TRUE){
	if(fold){
		nbin<-length(x)
		nchr<-(length(x)-1)*2
		tw<-sum(exp(x)[2:nbin])/sum(1/(1:(nchr-1)))
		return(tw)
		}
	if(!fold){
		nchr<-length(x)-1
		tw<-sum(exp(x)[2:nchr])/sum(1/(1:(nchr-1)))
		return(tw)
		}
	}
