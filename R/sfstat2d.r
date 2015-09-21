#simple stats from a 2d log-scaled sfs
#bp scales by pi and tw by number of sampled bases if SFS sums to 1 instead of bp

sfstat2d<-function(x,fold=FALSE,bp=3e7){

	x1<-log(rowSums(exp(x)))
	x2<-log(colSums(exp(x)))

	pi<-sfspi(x1,fold=fold)
	tw<-sfsw(x1,fold=fold)
	td<-tajimas(pi*bp,tw*bp,(length(x1)-1)*(fold+1))

	cat("row pi is ", pi, " ","tw is ", tw, " ", "tajima's D is ", td, "\n", sep="")

	pi<-sfspi(x2,fold=fold)
	tw<-sfsw(x2,fold=fold)
	td<-tajimas(pi*bp,tw*bp,(length(x2)-1)*(fold+1))

	cat("col pi is ", pi, " ","tw is ", tw, " ", "tajima's D is ", td, "\n", sep="")

	}

