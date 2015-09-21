#this awful function was meant to cover a few kinds of input files, but as long as col1=scaffold col2=position, you can specify a value to be plotted from any other column using 'col'. it will also plot a sliding window mean. there's lots of extraneous stuff in it that I was experimenting with. 

plot.fst<-function(dat, siglev=NULL, range=NULL, scaffold=NULL, type="fst", psize=.25, window=NULL, fai=NULL, col=3, mt="arithmetic"){
	
	if(!is.null(fai)){dat[,2]<-dat[,2]+index.lrt(dat, fai)}
	if(!is.null(scaffold)){dat[dat[,1]==scaffold,]->dat}
	if(!is.null(range)){dat[dat[,2]>=range[1]&dat[,2]<=range[2],]->dat}

	gpos<-"Genome position"
	
	if(type=="fst"){	
		yax<-c()
		if(col==5){yax<-"Fst"}
		if(col==4){yax<-"bg"; dat[,col]<-2*dat[,col]*(1-dat[,col])}
		if(col==3){yax<-"tg"; dat[,col]<-2*dat[,col]*(1-dat[,col])}
		}
	if(type=="ehh"){
		yax<-c()
		if(col==4){yax<-"iHS"}
		if(col==3){yax<-"xpEHH"}
		
		nans<-is.nan(dat[,3])|is.nan(dat[,4])
		dat<-dat[!nans,]
		
		}
	
	
	index.lrt<-function(dat, fai){
		offset<-as.numeric(fai[,3])
		names(offset)<-fai[,1]
		return(offset[dat[,1]])
		}

	if(type=="fst"){
		plot(x=dat[,2], y=dat[,col], pch=20, cex=psize, xlab=gpos, ylab=yax, ylim=c(0,1),main=scaffold)
		}else{plot(x=dat[,2], y=dat[,col], pch=20, cex=psize, xlab=gpos, ylab=yax, main=scaffold)}
	
	
	
	if(!is.null(siglev)){
		points(x=dat[,2][dat[,4]>(siglev)], y=dat[,col][dat[,4]>(siglev)], pch=20, cex=(psize*2), col="red")
		}
		
	if(!is.null(window)){
		if(mt=="arithmetic"){
			nwin<-(dat[length(dat[,1]),2]-dat[1,2])%/%window
			index<-dat[1,2]-1
			mid<-index+(.5*window)
			mat<-matrix(nrow=nwin, ncol=2, data=0)
			for(i in 1:nwin){
				mat[i,1]<-mid
				mat[i,2]<-mean(dat[dat[,2]>index&dat[,2]<=(index+window),col], na.rm=TRUE)
				index<-index+window
				mid<-mid+window
				}
			}
		if(mt=="harmonic"){
			nwin<-(dat[length(dat[,1]),2]-dat[1,2])%/%window
			index<-dat[1,2]-1
			mid<-index+(.5*window)
			mat<-matrix(nrow=nwin, ncol=2, data=0)
			for(i in 1:nwin){
				mat[i,1]<-mid
				mat[i,2]<-1/mean(1/(dat[dat[,2]>index&dat[,2]<=(index+window),col]), na.rm=TRUE)
				index<-index+window
				mid<-mid+window
				}
			}

			lines(x=mat[,1], y=mat[,2], col="green", lwd=2)
		}
		
	}


