pullMLgenos<-function(scaf,vcf="~/popgen/variants/bowfree/ALL1/hallsnps.vcf.gz"){

	cline<-paste("tabix ", vcf, " ", scaf,sep="")
	header<-scan(pipe(paste("tabix -H ",vcf," | tail -n 1 ",sep=" ")),what="character")
	gts<-read.table(pipe(cline),stringsAsFactors=FALSE)

	colnames(gts)<-header
	gts2<-as.matrix(gts[,10:393])
	gts2<-gsub(".*:","",gts2)

	mlgt<-function(y){
		if(y=="."){return(NA)}
		(which(unlist(strsplit(x=y,split=","))==0)-1)[1]
		}

	vmlgt<-function(z){
		sapply(z,FUN=mlgt,USE.NAMES=FALSE)
		}

	for(i in 1:length(gts2[1,])){
		tmp<-vmlgt(gts2[,i])
		gts2[,i]<-tmp

		}

	gts2<-cbind(gts[,2],gts2)
	class(gts2)<-"numeric"
	return(gts2)
	}

