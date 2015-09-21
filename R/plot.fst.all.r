#this function should plot Fst across the whole genome, coloring scaffolds alternately "black" and "darkgray". it subsamples fst values every "sub" in order to prevent the plot from becoming too insane.


plot.fst.all<-function(dat,sub=100){

	bgcol<-c("black","darkgray")
	subsamp<-seq(1,length(dat[,1]),sub)
	dat<-dat[subsamp,]

	datcol<-c()
	dattab<-cbind(table(dat[,1]),bgcol)[unique(dat[,1]),]
	for(i in 1:length(dattab[,1])){
		datcol<-c(datcol, rep(dattab[i,2],dattab[i,1]))
		}

	plot(dat[,5],col=datcol,pch=20,cex=.1,ylim=c(0,1))

	}
