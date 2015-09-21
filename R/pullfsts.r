#this function pulls Fst values from a given scaffold for a bgzipped and tabix-indexed table containing Fst for all snps in the genome. 
	# in this case it pulls the values for four population pairs and puts each in a list. 
	
pullfsts<-function(scaffold, fst1="BI.NBH",fst2="BP.F",fst3="NYC.SH",fst4="ER.KC"){

	res<-list()
	res[[fst1]]<-read.table(pipe(paste("tabix /home/nreid/popgen/variants/bowfree/ALL1/", fst1,".wcfst.gz ", scaffold, collapse="", sep="")))
	res[[fst2]]<-read.table(pipe(paste("tabix /home/nreid/popgen/variants/bowfree/ALL1/", fst2,".wcfst.gz ", scaffold, collapse="", sep="")))
	res[[fst3]]<-read.table(pipe(paste("tabix /home/nreid/popgen/variants/bowfree/ALL1/", fst3,".wcfst.gz ", scaffold, collapse="", sep="")))
	res[[fst4]]<-read.table(pipe(paste("tabix /home/nreid/popgen/variants/bowfree/ALL1/", fst4,".wcfst.gz ", scaffold, collapse="", sep="")))

	return(res)
	}

