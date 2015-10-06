#feed this function a set of windows and corresponding statistics passing a threshold. 
	#it will merge them if they overlap or are within 'buffer' of each other and calculate stats. 
	#modified from a previous version in process_stats.R
	#
mergewin<-function(win, stat, qu, tails=TRUE, buff=5000){
	
	stat<-as.numeric(stat)
	stat[is.na(stat)]<-0
	if(tails=="lesser"){
		qu<-qu[1]
		stat<-stat-qu
		stat[stat>0]<-0	
		}
	if(tails=="greater"){
		qu<-qu[length(qu)]
		stat<-stat-qu
		stat[stat<0]<-0
		}
	if(tails=="both"){
		if(length(qu)==1){stop("two quantiles must be provided for option 'both'")}
		stat[stat>qu[1]&stat<qu[2]]<-0
		stat[stat<qu[1]]<-stat[stat<qu[1]]-qu[1]
		stat[stat>qu[2]]<-stat[stat>qu[2]]-qu[2]
		}
	
	#chrom,start,end,length,auc,highpoint,count
	out<-cbind(win[1,],win[1,3]-win[1,2],stat[1],stat[1],1)
	names(out)<-c("chrom","start","end","length","AUC","highpoint","count")
	for(i in 2:length(win[,1])){
		
		ind<-length(out[,1])
		if(out[ind,1]==win[i,1]&(out[ind,3]+buff)>=win[i,2]){
			out[ind,3]<-win[i,3]
			out[ind,4]<-out[ind,3]-out[ind,2]
			out[ind,5]<-out[ind,5]+stat[i]
			if(abs(out[ind,6])<abs(stat[i])){out[ind,6]<-stat[i]}
			if(stat[i]!=0){out[ind,7]<-out[ind,7]+1}
			}else{
					out<-rbind(out,setNames(cbind(win[i,],win[i,3]-win[i,2],stat[i],stat[i],1),names(out)))
					}		
		}
	return(out)
	}
	
	