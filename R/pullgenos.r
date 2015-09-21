#this function pulls genotypes from a tabix-indexed vcf containing biallelic snps
#it tosses most of the information and produces a table of positions, and genotypes coded as 0,1,2

pullgenos<-function(scaf,vcf="~/popgen/variants/bowfree/ALL1/hallsnps.vcf.gz"){

	cline<-paste("tabix ", vcf, " ", scaf,sep="")
	header<-scan(pipe(paste("tabix -H ",vcf," | tail -n 1 ",sep=" ")),what="character")
	gts<-read.table(pipe(cline),stringsAsFactors=FALSE)

	colnames(gts)<-header
	gts2<-as.matrix(gts[,10:393])
	gts2<-gsub(":.*","",gts2)
	gts2[gts2=="."]<-NA
	gts2[gts2=="0/0"]<-0
	gts2[gts2=="0/1"]<-1
	gts2[gts2=="1/1"]<-2

	gts2<-cbind(gts[,2],gts2)
	class(gts2)<-"numeric"
	return(gts2)
	}

