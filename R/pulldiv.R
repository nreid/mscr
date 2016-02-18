pulldiv <- function(scaf,stat="pi",win="5kb1kbwindow",popids=c("bi","bp","er","f","kc","nbh","nyc","sh"),path="~/popgen/variants/bowfree/angsd/"){
	pops <- popids
	#tabix sh.pi.5kb1kbwindow.test.gz Scaffold9893
	div <- list()
	for(i in 1:length(pops)){

		fn <- paste(pops[i],stat,win,"test.gz",sep=".",collapse=".")
		print(fn)
		fn <- paste(path,fn,sep="",collapse="")
		print(fn)
		cline <- paste("tabix ",fn,scaf,sep=" ",collapse=" ")
		print(cline)
		div[[i]] <- read.table(pipe(cline),stringsAsFactors=FALSE)
		
		}

	out <- div[[1]][,c(2:3,6)]
	for(i in 1:length(pops)){

		out <- cbind(out,as.numeric(as.character(div[[i]][,4])))

		}

	colnames(out) <- c("start","end","nsites",pops)

	return(out)
	}