##reads a region of a tabix indexed vcf (using tabix)
##outputs a table with VCF colnames

read.vcf <- function(vcf="~/popgen/variants/bowfree/ALL1/hallsnps.vcf.gz",reg){

  cl1 <- paste("tabix",vcf,reg,sep=" ")
  cl2 <- paste("zcat", vcf," | head -n 1000 | grep '^.CHROM' ")
  con <- pipe(cl1)
  scaf <- read.table(con,stringsAsFactors=FALSE)
  con <- pipe(cl2)
  hed <- scan(con,what="character")
  close(con)
  colnames(scaf)<-hed
  return(scaf)
  }
