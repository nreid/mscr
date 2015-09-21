##this function parses output from a single population in ms.
#it outputs a table of summary statistics

msparse<-function(x,nsites){

	he<-grep("segsites", x)
	ind<-he[2]-he[1]-3

	out<-vector("list", ndat)
	out2<-vector("list", ndat)

 	piest<-function(x){

		inds<-dim(x)[1]
		sum(colSums(x)*(inds-colSums(x))/choose(inds,2))

		}

	twest<-function(x){

		inds<-dim(x)[1]
		sum(colSums(x)>0)/sum(1/(1:(inds-1)))

		}

	for(i in 1:ndat){

		if(x[he[i]]=="segsites: 0"){
			out[[i]]<-0
			out2[[i]]<-0
			}
			else{

				outmat<-x[(he[i]+2):(he[i]+1+ind)]
				outmat<-do.call(rbind,strsplit(outmat,""))
				class(outmat)	<-"numeric"
				out2[[i]]<-twest(outmat)
				out[[i]]<-piest(outmat)
					}
		}

	cbind(unlist(out),unlist(out2))

	}


