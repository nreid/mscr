
#this function pulls genotypes from a tabix-indexed vcf containing biallelic snps that have been phased
#it tosses most of the information and produces a table of positions, and haplotypes coded as 0,1

pullgenos.phased <- function(scaf,vcf="~/popgen/variants/bowfree/ALL1/hallsnps.imp.vcf.gz"){

	cline<-paste("tabix ", vcf, " ", scaf,sep="")
	header<-scan(pipe(paste("tabix -H ",vcf," | tail -n 1 ",sep=" ")),what="character")
	gts<-read.table(pipe(cline),stringsAsFactors=FALSE)

	colnames(gts)<-header
	gts2<-as.matrix(gts[,-c(1:9)])
	gts2<-gsub(":.*","",gts2)
	gts3 <- c()
	for(i in 1:dim(gts2)[2]){
		temp <- str_split(gts2[,i],"\\|")
		temp <- do.call(rbind,temp)
		colnames(temp) <- paste(colnames(gts2)[i],c("_a","_b"),sep="")
		gts3 <- cbind(gts3,temp)
		}
	gts3<-cbind(gts[,2],gts3)
	class(gts3)<-"numeric"
	return(gts3)

	}