#pi, tw, tajimas from a log-scaled sfs

sfstat<-function(x,fold=FALSE,bp=3e7){
	pi<-sfspi(x,fold=fold)
	tw<-sfsw(x,fold=fold)
	td<-tajimas(pi*bp,tw*bp,(length(x)-1)*(fold+1))

	cat("pi is ", pi, "\n","tw is ", tw, "\n", "tajima's D is ", td, "\n", sep="")
	}
