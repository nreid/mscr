pimage<-function(x,pop=NULL,slo=NULL,col=c("yellow","green","blue")){
  x<-x[,-1]
  if(!is.null(spo)){
    spo<-grep(pop,colnames(x))
    x<-x[,spo]
  }
  if(!is.null(slo)){x<-x[slo,]}

  image(x,col=col)

  }
