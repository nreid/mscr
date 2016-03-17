#this function pulls Fst values from a given scaffold for a bgzipped and tabix-indexed table containing Fst for all snps in the genome. 
	# in this case it pulls the values for four population pairs and puts each in a list. 
	
pullfsts_win<-function(scaffold, fst1="BI.NBH",fst2="BP.F",fst3="NYC.SH",fst4="ER.KC"){

	res<-list()
	res[[fst1]]<-read.table(pipe(paste("tabix /home/nreid/popgen/variants/bowfree/angsd_fst/", fst1,".fst.win.gz ", scaffold, collapse="", sep="")))
	res[[fst2]]<-read.table(pipe(paste("tabix /home/nreid/popgen/variants/bowfree/angsd_fst/", fst2,".fst.win.gz ", scaffold, collapse="", sep="")))
	res[[fst3]]<-read.table(pipe(paste("tabix /home/nreid/popgen/variants/bowfree/angsd_fst/", fst3,".fst.win.gz ", scaffold, collapse="", sep="")))
	res[[fst4]]<-read.table(pipe(paste("tabix /home/nreid/popgen/variants/bowfree/angsd_fst/", fst4,".fst.win.gz ", scaffold, collapse="", sep="")))
	for(i in 1:4){res[[i]][,4] <- as.numeric(as.character(res[[i]][,4]))}
	for(i in 1:4){res[[i]][,5] <- as.numeric(as.character(res[[i]][,5]))}
	for(i in 1:4){res[[i]][,6] <- as.numeric(as.character(res[[i]][,6]))}
	return(res)
	}

