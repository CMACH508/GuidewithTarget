args=commandArgs(T)
file1 = args[1] # POI_Target_Pair:ago3_POI_Target_Pair.txt
data_add = args[2]
reference = args[3]
iteration = args[4]
cellline = unlist(strsplit(file1,"_"))[1]
reference = ifelse(reference=="1","transcriptome","genome")
options(scipen=999)

library(dplyr)
library("data.table")
library(parallel)




#file1 = "ago3_POI_Target_Pair.txt"
#data_add = "/data3/zenglin/temp/LinZeng/degradome/DrosoDataStore/ago3/"
#reference = "2"
#cellline = "ago3"
#reference = ifelse(reference=="1","transcriptome","genome")
result <- readLines(paste0(data_add,file1))
rownum_query <- which(startsWith(result,"Query="))


#system(paste0("rm ",data_add,cellline,"_Parsed_NPOI_",iteration,".txt"))
#system(paste0("rm ", data_add,cellline,"_Parsed_Filter_NPOI_",iteration,".txt"))

ParseEachQuery <- function(query_num){

  removesym <- function(i){
  	pro_i <- substr(i,3,nchar(i))
  	  return(pro_i)
  }
  removelast <- function(i){
  	pro_i <- substr(i,1,nchar(i)-3)
  	return(pro_i)
  }
  extractpos <- function(i){
    b1 <- as.numeric((unlist(strsplit(i,"::"))[1] %>% strsplit("\\|") %>% unlist)[3])
    return(b1)
  }
  extrastart <- function(i){
    s1 <- unlist((unlist(strsplit(i,"::"))[2] %>% strsplit(":") %>% unlist)[2] %>% strsplit("-"))[1]
    return(as.numeric(s1))
  }
  df2 <- data.frame(query_id=character(), sub_id=character(), q_start=character(),
    q_end=character(), s_start=character(),  s_end=character(),
    q_match_pos=character(), q_mismatch_pos=character(),
    identity=character(), alignment_length = numeric(),  match_num = numeric(),
    mismatch_num = numeric(), stringsAsFactors=FALSE)
    
  query_row = rownum_query[query_num]
  query_extrat <- result[query_row:(rownum_query[query_num+1]-1)]
  sub_id_list = query_extrat[which(startsWith(query_extrat,"> "))]
  query_id = substr(query_extrat[1],8,nchar(query_extrat[1]))
  Nohist = query_extrat[6]
  
	if(Nohist=="***** No hits found *****"){
		print("No Hits")
	}else{
    for(subject in 1:length(sub_id_list)){
      sub_id = sub_id_list[subject]
      sub_rownum = which(startsWith(query_extrat,sub_id))
		  sub_rownum1 = sub_rownum+11
      content_rel <- query_extrat[sub_rownum:sub_rownum1]
      identity <- (unlist(strsplit(content_rel[5],"), Gaps = "))[1] %>% strsplit("\\(") %>% unlist)[2]
      alignment_summary = (unlist(strsplit(content_rel[5]," Identities = "))[2] %>% strsplit(" \\(") %>% unlist)[1]
      alignment_length = as.numeric(unlist(strsplit(alignment_summary,"/"))[2])
      match_num = as.numeric(unlist(strsplit(alignment_summary,"/"))[1])
      mismatch_num = alignment_length - match_num
      temp1 = unlist(strsplit(content_rel[8]," "))
      temp1 = temp1[which(temp1!="")]
      q_start <- temp1[2]
		  q_end <- temp1[4]
		  temp2 = unlist(strsplit(content_rel[10]," "))
		  temp2 = temp2[which(temp2!="")]
		  s_end <- temp2[2]
		  s_start <- temp2[4]
		  matchpos <- content_rel[9] #"           ||| || ||||||||||||||"
		  q_match_pos = ""
		  q_mismatch_pos = ""
		  count_padding = 0
      while(substr(matchpos,count_padding+1,count_padding+1)!="|"){
        		count_padding = count_padding + 1
		  }
	  	for(i in 1:nchar(matchpos)){
			      #print("==============")
        		symbol <- substr(matchpos,i,i)
       			if(symbol=="|"){
		                q_match_pos = c(q_match_pos,i-count_padding+as.numeric(q_start)-1)
        		}
		  }
      q_match_pos = q_match_pos[2:length(q_match_pos)]
	    q_mismatch_pos = setdiff(c(as.numeric(q_start):(as.numeric(q_end)-1)),q_match_pos)
	    q_match_pos <- paste0(q_match_pos,collapse=",")
    	q_mismatch_pos <- paste0(q_mismatch_pos,collapse=",")
	    pairresult <- data.frame(query_id=query_id, sub_id=sub_id,q_start=q_start,q_end=q_end,s_tart=s_start,s_end=s_end,q_match_pos=q_match_pos,q_mismatch_pos=q_mismatch_pos,alignment_length=alignment_length,
      match_num=match_num,mismatch_num=mismatch_num,identity=identity)
	    df2 <- rbind(df2,pairresult)
    }
  }
  df2$sub_id <- unlist(lapply(df2$sub_id,removesym))
  df2$dis1 <- 10-as.numeric(df2$q_start)
  df2$pos1_based <- unlist(lapply(df2$sub_id,extractpos))
  df2$start_abs <- unlist(lapply(df2$sub_id,extrastart))
  relpos <- df2$pos1_based - df2$start_abs
  dis <- as.numeric(df2$s_end) - relpos
  df2$dis <- dis
  df2$relpos <- relpos
  truedf2 <- df2[which(df2$dis==df2$dis1),] 
  truedf2 <- truedf2[which(truedf2$dis > 0),]
  write.table(file=paste0(data_add,cellline,"_Parsed_NPOI_",iteration,".txt"),df2,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,sep="@")
  write.table(file=paste0(data_add,"blast_result/",cellline,"_Parsed_Filter_NPOI_",query_num,".txt"),truedf2,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE,sep="@")
}


no_of_cores <- 2
cl <- makePSOCKcluster(no_of_cores)
clusterExport(cl, c("rownum_query","result","data_add","cellline","iteration"))
clusterEvalQ(cl,library("dplyr"))
results <- parLapply(cl,1:(length(rownum_query)-1),ParseEachQuery)
stopCluster(cl)
#gc()
