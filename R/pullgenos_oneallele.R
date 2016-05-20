#this function pulls genotypes from a tabix-indexed vcf containing biallelic snps
#it tosses most of the information and produces a table of positions, and genotypes coded as 0,1,2

pullgenos_oneallele<-function(scaf,vcf="~/popgen/variants/bowfree/ALL1/hallsnps.vcf.gz"){

	cline<-paste("tabix ", vcf, " ", scaf,sep="")
	header<-scan(pipe(paste("tabix -H ",vcf," | tail -n 1 ",sep=" ")),what="character")
	gts<-read.table(pipe(cline),stringsAsFactors=FALSE)

	sam<-function(x){
		x <- unlist(str_split(x,":"))[c(3,5)]
		x <- as.numeric(x)!=0
		if(any(is.na(x))){return(NA)}
		if(sum(x)==0|sum(x)==2){return(NA)}
		x <- x/sum(x)
		x <- sample(x=c(0,1),size=1,prob=x)
		return(x)
		}
	samcol <- function(x){
		sapply(x,FUN=sam)
		}


	colnames(gts)<-header
	gts2<-as.matrix(gts[,10:393])

	gts2 <- apply(gts2,MAR=2,FUN=samcol)

	gts2<-cbind(gts[,2],gts2)
	class(gts2)<-"numeric"
	return(gts2)
	}

