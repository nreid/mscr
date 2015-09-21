#this function parses ms output generated from two populations.
#it outputs a table of summary statistics
#it requires the function 'wc' from hierfstat to caculate fst
#n1 and n2 are numbers of 'chromosomes' simulated
#x is output from ms read in using scan(x,what="character",sep="\n")

msparse.2pops<-function(x,n1,n2,nsites){

	he<-grep("segsites", x)
	ndat<-length(he)
	ind<-he[2]-he[1]-3

	out1pi<-vector("list", ndat)
	out1tw<-vector("list", ndat)
	out1tj<-vector("list", ndat)
	out2pi<-vector("list", ndat)
	out2tw<-vector("list", ndat)
	out2tj<-vector("list", ndat)
	outfst<-vector("list", ndat)


 	piest<-function(x){

		inds<-dim(x)[1]
		sum(colSums(x)*(inds-colSums(x))/choose(inds,2))

		}

	twest<-function(x){

		inds<-dim(x)[1]
		sum(colSums(x)>0)/sum(1/(1:(inds-1)))

		}

	fstest<-function(x,num1,num2){
		x<-x+10
		pops<-c(rep(1,num1),rep(2,num2))
		x<-cbind(pops,x)
		wc(x,diploid=FALSE)$FST
		}

	for(i in 1:ndat){

		if(x[he[i]]=="segsites: 0"){
			out1pi[[i]]<-0
			out1tw[[i]]<-0
			out1tj[[i]]<-0
			out2pi[[i]]<-0
			out2tw[[i]]<-0
			out2tj[[i]]<-0
			outfst[[i]]<-0

			}
			else{

				outmat<-x[(he[i]+2):(he[i]+1+ind)]
				outmat<-do.call(rbind,strsplit(outmat,""))
				class(outmat)	<-"numeric"
				out1pi[[i]]<-piest(outmat[1:n1,])
				out1tw[[i]]<-twest(outmat[1:n1,])
				out1tj[[i]]<-tajimas(out1pi[[i]],out1tw[[i]],n1)
				out2pi[[i]]<-piest(outmat[(n1+1):(n1+n2),])
				out2tw[[i]]<-twest(outmat[(n1+1):(n1+n2),])
				out2tj[[i]]<-tajimas(out2pi[[i]],out2tw[[i]],n2)
				outfst[[i]]<-fstest(outmat,n1,n2)
					}
		if(ndat%%10==0){print(i)}
		}

	cbind(unlist(out1pi),unlist(out1tw),unlist(out1tj),unlist(out2pi),unlist(out2tw),unlist(out2tj),unlist(outfst))

	}

