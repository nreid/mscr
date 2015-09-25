###see file Frequency_LD_scripts.txt for some tests

###functions to estimate allele frequencies given precalculated genotype likelihoods from Lynch 2008

#Lpgl gives the likelihood of an allele frequency, p, given a vector of genotype likelihoods, gl
#x is a vector of 3 genotype likelihoods
#p is an allele frequency
#assume p is alt allele, and gls are ordered (1-p)^2,2*(1-p)*p,p^2
#gls are from freebayes, log10 scaled. 



Lpgl_ind<-function(gl,p){
	
	if(sum(is.na(gl))){
		lik<-1
		}
	else{
		gl<-10^gl
		lik<-((p^2)*gl[3])+(2*p*(1-p)*gl[2])+(((1-p)^2)*gl[1])
		}
	lik<-log(lik)		
	}
	

#function to apply function Lpgl across a table of genotype likelihoods for multiple individuals at a single site
	#columns are gls, ordered (1-p)^2,2*(1-p)*p,p^2, rows are individuals
	
Lpgl_site<-function(pt,glt){
	
	lik<-apply(glt,MAR=1,FUN=Lpgl_ind,p=pt)
	lik<-sum(lik)
	return(lik)
	
	}
	
#function to calculate frequencies 
#vcf is vcf file read in using read.table()
#sub is a vector of individual numbers (not column numbers) to calculate the frequency from

calcfreqs<-function(vcf,sub=NULL){
	
	nsites<-length(vcf[,1])
	if(is.null(sub)){
		sub<-1:(length(vcf[1,])-9)
		}
		
	vcf<-as.matrix(vcf[,sub+9])
	vcf<-gsub(".*:","",vcf)
	
	
	freqs<-rep(NA,nsites)
	for(i in 1:nsites){
		
		gltab<-do.call(rbind,strsplit(x=vcf[i,],split=","))
		gltab[gltab=="."]<-NA
		class(gltab)<-"numeric"
		freqs[i]<-optimize(f=Lpgl_site, interval=c(0,1),glt=gltab,maximum=TRUE)$maximum
		if(i%%1000==0){
			cat(i, " iterations are done\n")
			}		
		}
	
	return(freqs)
	
	}


###this function takes genotype likelihoods from two loci gl1,gl2
###a value of D, and allele frequencies at each locus, p,q
###and returns the natural log scaled likelihood of D. 

Ld_glpq_ind<-function(gl1,gl2,d,p,q){
	
	if(any(is.na(c(gl1,gl2)))){return(log(1))}
	altalt<-log(p*q+d,base=10)
	altref<-log(p*(1-q)-d,base=10)
	refalt<-log((1-p)*q-d,base=10)
	refref<-log((1-p)*(1-p)+d,base=10)
	
	lscores<-c(
		gl1[3]+gl2[3]+altalt+altalt,
		gl1[3]+gl2[1]+altref+altref,
		gl1[1]+gl2[3]+refalt+refalt,
		gl1[1]+gl2[1]+refref+refref,
		gl1[3]+gl2[2]+altalt+altref+log(2,base=10),	
		gl1[1]+gl2[2]+refalt+refref+log(2,base=10),	
		gl1[2]+gl2[3]+altalt+refalt+log(2,base=10),	
		gl1[2]+gl2[1]+altref+refref+log(2,base=10),
		gl1[2]+gl2[2]+altalt+refref+log(2,base=10),	
		gl1[2]+gl2[2]+altref+refalt+log(2,base=10)
		)	
	log(sum(10^lscores))
	}


##calculate likelihood of D given gls and frequencies for two sites for all individuals
Ld_glpq_sitepair<-function(ds,gl1s,gl2s,ps,qs){
	
	
	l_ind<-apply(cbind(gl1s,gl2s),MAR=1,FUN=function(x){Ld_glpq_ind(x[1:3],x[4:6],d=ds,p=ps,q=qs)})
	
	sum(l_ind)
	
	}



calcD<-function(vcf,sub=NULL,pvec=NULL,maxdist=NULL,bypos=FALSE){
	
	if(is.null(pvec)){
		print("calculating allele frequencies")
		pvec<-calcfreqs(vcf=vcf,sub=sub)
		}
	subsites<-pvec>.2&pvec<.8
	pvec<-pvec[subsites]
	vcf<-vcf[subsites,]
	
	pos<-vcf[,2]
	nsites<-dim(vcf)[1]

	if(is.null(sub)){
		sub<-1:(length(vcf[1,])-9)
		}		
	vcf<-as.matrix(vcf[,sub+9])
	vcf<-gsub(".*:","",vcf)
	
	
	allpairs<-t(combn(1:nsites,2))
	if(!is.null(maxdist&!bypos)){
		allpairs<-allpairs[allpairs[,2]-allpairs[,1]<maxdist,]
		}	
	if(!is.null(maxdist&bypos)){
		pairpos<-t(combn(pos,2))
		subpos<-pairpos[,2]-pairpos[,1]<maxdist
		allpairs<-allpairs[subpos,]
		}	


	npairs<-dim(allpairs)[1]
	out<-matrix(nrow=npairs,ncol=5)
	
	cat("calculating D values\n")
	for(i in 1:npairs){
		
		GL1<-do.call(rbind,strsplit(split=",",vcf[allpairs[i,1],]))
		GL2<-do.call(rbind,strsplit(split=",",vcf[allpairs[i,2],]))
		class(GL1)<-"numeric"
		class(GL2)<-"numeric"
		P1<-pvec[allpairs[i,1]]
		P2<-pvec[allpairs[i,2]]
		dmin<-(-min(P1*P2,(1-P1)*(1-P2)))
		dmax<-min(P1*(1-P2),(1-P1)*(P2))
		out[i,1]<-pos[allpairs[i,1]]
		out[i,2]<-pos[allpairs[i,2]]
		out[i,3]<-P1
		out[i,4]<-P2
		out[i,5]<-optimize(f=Ld_glpq_sitepair, interval=c(dmin,dmax),gl1s=GL1,gl2s=GL2,ps=P1,qs=P2,maximum=TRUE)$maximum
		
		if(i%%1000==0){
			cat(i, " iterations are done\n")
			}		

		}
	
	return(out)
	
	}



