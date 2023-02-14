args=commandArgs(T) 
file1 = args[1] #  "sptd_total_pi6.txt"
data_add = args[2] #  data_add
cell_line = args[3] # sptd
data_end = as.numeric(args[4])
reference = args[5]
reference = ifelse(reference=="1","Transcriptome","Genome")
options(scipen = 100)

#file_pi6 <- "sptd_total_pi6.txt"
#data_add = "/home/user/Zenglin/degradome/datastore_final/sptd/"
#setwd(data_add)
#cell_line = "sptd"
"%notin%" <- Negate("%in%")
degwt_extract <- function(i){
  i_proc <- as.numeric(unlist(strsplit(i,"\\|"))[2])
  return(i_proc)
}


filelist <- dir(data_add)

filelist <- filelist[startsWith(filelist,paste0(cell_line,"_Parsed_Filter_NPOI_"))]
#print(filelist)

df_sum <- data.frame(grp=character(),
                     frac=numeric(),
                     start=numeric(),
                     end=numeric(),
                     cat=character(),
                     X=integer(),
                     iteration=numeric(),
                     stringsAsFactors=FALSE)

for(iter in 1:length(filelist)){
  
  nonpi6 <- filelist[iter]
  print(nonpi6)
  df_result <- data.frame(X=integer(),
                          grp=character(),
                          fract=numeric(),
                          pattern=character(),
                          cat_start=character(),
                          cat_end=character(),
                          stringsAsFactors=FALSE)
  
  total_pi6 <- read.table(paste0(data_add,file1),sep="@",header=F)
  total_nonpi6 <- read.table(paste0(data_add,nonpi6),sep="@",header=F)
  name2sign <- c("query_id","sub_id","q_start","q_end","s_tart","s_end",
                 "q_match_pos","q_mismatch_pos","alignment_length","match_num",
                 "mismatch_num","identity","dis1","pos1_based","start_abs","dis","relpos")
  colnames(total_pi6) = name2sign
  colnames(total_nonpi6) = name2sign
  
  total_pi6$deg_wt <- unlist(lapply(total_pi6$sub_id,degwt_extract))
  total_nonpi6$deg_wt <- unlist(lapply(total_nonpi6$sub_id,degwt_extract))
  
  pi6fullmatch <- total_pi6[which(total_pi6$identity=="100%"),]
  pi6fullmatch <- pi6fullmatch[which(pi6fullmatch$q_mismatch_pos==""),]
  pi6fullmatch$grp = NA
  pi6fullmatch$frac = 0 
  pi6fullmatch$q_match_pos = gsub("^1,","",pi6fullmatch$q_match_pos)
  #pi6fullmatch$deg_wt = 0
  for(rownum in 1:nrow(pi6fullmatch)){
    q_start <- unlist(strsplit(pi6fullmatch[rownum,"q_match_pos"],","))[1]
    q_end <- pi6fullmatch[rownum,"q_end"]
    q_end <- ifelse(q_end>data_end, data_end, q_end)
    grp <- paste0(q_start,"_",q_end)
    pi6fullmatch[rownum,"grp"] = grp
    deg_wt <- as.numeric(unlist(strsplit(pi6fullmatch[rownum,"sub_id"],"\\|"))[2])
    #pi6fullmatch[rownum,"deg_wt"] = deg_wt
  }
  require(dplyr)
  pi6stat <- pi6fullmatch %>%
    select(deg_wt, grp) %>% 
    group_by(grp) %>%
    summarize(frac = sum(deg_wt))
  temp <- c()
  for(i in 2:9){
    for(j in 9:data_end){
      tem = paste0(i,"_",j)
      temp <- c(temp,tem)
    }
  }
  temp <- temp[which(temp %notin% pi6stat$grp)]
  comp_stat <- data.frame(grp=temp,frac=rep(0,length(temp)))
  pi6stat <- rbind(pi6stat,comp_stat)
  pi6stat <- pi6stat[order(pi6stat$grp),]
  pi6stat$start <- as.numeric(unlist(strsplit(pi6stat$grp,"_"))[seq(1,nrow(pi6stat)*2,2)])
  pi6stat$end <- as.numeric(unlist(strsplit(pi6stat$grp,"_"))[seq(2,nrow(pi6stat)*2,2)])
  pi6stat <- pi6stat[order(pi6stat$start,pi6stat$end),]
  pi6stat$frac <- pi6stat$frac/(sum(total_pi6$deg_wt)+sum(total_nonpi6$deg_wt))
  
  #===================================
  
  nonpi6fullmatch <- total_nonpi6[which(total_nonpi6$identity=="100%"),]
  #nonpi6fullmatch <- nonpi6fullmatch[which(nonpi6fullmatch$q_mismatch_pos==""),]
  nonpi6fullmatch$grp = NA
  nonpi6fullmatch$frac = 0
  nonpi6fullmatch$q_match_pos = gsub("^1,","",nonpi6fullmatch$q_match_pos)
  for(rownum in 1:nrow(nonpi6fullmatch)){
    q_start <- unlist(strsplit(nonpi6fullmatch[rownum,"q_match_pos"],","))[1]
    q_end <- nonpi6fullmatch[rownum,"q_end"]
    q_end <- ifelse(q_end>data_end, data_end, q_end)
    grp <- paste0(q_start,"_",q_end)
    nonpi6fullmatch[rownum,"grp"] = grp
  }
  require(dplyr)
  nonpi6stat <- nonpi6fullmatch %>%
    select(deg_wt, grp) %>% 
    group_by(grp) %>%
    summarize(frac = sum(deg_wt))
  temp <- c()
  for(i in 2:9){
    for(j in 9:data_end){
      tem = paste0(i,"_",j)
      temp <- c(temp,tem)
    }
  }
  temp <- temp[which(temp %notin% nonpi6stat$grp)]
  comp_stat <- data.frame(grp=temp,frac=rep(0,length(temp)))
  nonpi6stat <- rbind(nonpi6stat,comp_stat)
  nonpi6stat <- nonpi6stat[order(nonpi6stat$grp),]
  nonpi6stat$start <- as.numeric(unlist(strsplit(nonpi6stat$grp,"_"))[seq(1,nrow(nonpi6stat)*2,2)])
  nonpi6stat$end <- as.numeric(unlist(strsplit(nonpi6stat$grp,"_"))[seq(2,nrow(nonpi6stat)*2,2)])
  nonpi6stat <- nonpi6stat[order(nonpi6stat$start,nonpi6stat$end),]
  nonpi6stat$frac <- nonpi6stat$frac/(sum(total_pi6$deg_wt)+sum(total_nonpi6$deg_wt))
  
  combi <- data.frame(grp=pi6stat$grp,frac=pi6stat$frac - nonpi6stat$frac,
                      start=pi6stat$start, end=pi6stat$end)
  
  combi$cat <- "combi"
  nonpi6stat$cat <- "nonpi6"
  pi6stat$cat <- "pi6"
  result <- rbind(combi,nonpi6stat,pi6stat)
  result$X <- rep(c(1:length(rep(c(9:data_end),9-2+1))),3)
  print(iter)
  result$iteration = iter
  
  df_sum <- rbind(result,df_sum)
}


df_sum$iteration = as.factor(df_sum$iteration)

# p <- ggplot(df_sum[which(df_sum$start==2 & df_sum$cat=="combi"),], 
#             aes(x=end, y=frac, group = iteration, color=iteration))+ 
#   geom_boxplot()+
#   #geom_line(aes(linetype=cat)) + 
#   #geom_point(aes(shape=iteration))+
#   theme_classic()+
#   #geom_vline(xintercept=c((data_end-9+1),(data_end-9+1)*2,(data_end-9+1)*3,(data_end-9+1)*4,(data_end-9+1)*5,(data_end-9+1)*6,(data_end-9+1)*7,(data_end-9+1)*8),linetype='dashed')+
#   theme(plot.title = element_text(hjust=0.5, color="black", size=20, face="bold"))+
#   scale_color_brewer(palette="Reds")

library(ggplot2)
p <- ggplot(df_sum[which(df_sum$cat=="combi"),], aes(x=end, y=frac, group = iteration, color=iteration))+ 
  #geom_line(aes(linetype=cat)) + 
  geom_point(size=1)+
  #geom_boxplot(df_sum[which(df_sum$start==9 & df_sum$cat=="combi"),], aes(x=end, y=frac, group = iteration, color=iteration))+
  theme_bw()+
  #geom_vline(xintercept=c((data_end-9+1),(data_end-9+1)*2,(data_end-9+1)*3,(data_end-9+1)*4,(data_end-9+1)*5,(data_end-9+1)*6,(data_end-9+1)*7,(data_end-9+1)*8),linetype='dashed')+
  theme(plot.title = element_text(hjust=0.5, color="black", size=20, face="bold"),legend.position = "none")+
  #scale_color_brewer(palette="Reds")+
  facet_grid(start ~ .)+
  scale_x_continuous(breaks = c(9:data_end))

ggsave(paste0(data_add,"FullMatch_Scatter_",length(filelist),"_",data_end,"_",reference,".pdf"),dpi=300,width = 6,height = 10)
write.table(file = paste0(data_add,"FullMatch_Scatter_",length(filelist),"_",data_end,"_",reference,".txt"),
            df_sum,quote=FALSE,row.names=FALSE,sep="\t")

mat = as.data.frame(matrix(0,(data_end-9+1), (9-2+1)))
rownames(mat) = c(9:data_end)
colnames(mat) = c(2:9)
for(row_num in rownames(mat)){
  for(col_num in colnames(mat)){
    df_temp <- df_sum[which(df_sum$start==as.numeric(col_num) & 
                              df_sum$end==as.numeric((row_num))&
                              df_sum$cat=="combi"),]
    value_temp <- median(df_temp$frac)
    mat[row_num,col_num] = value_temp
  }
}


colnames(mat) <- c(paste0("g",c(2:9),"-gX"))
library("ComplexHeatmap")
library("circlize")
plt2 <- Heatmap(as.matrix(mat), 
                name = "fcleav",
                km=1,
                #clustering_distance_columns = distype[1],
                #clustering_method_columns = methtype[1],
                colorRamp2(c(min(mat),max(mat)), c("white","red")),
                show_column_names = T,
                cluster_columns = F,
                cluster_rows = F,
                show_row_names = T,
                row_names_gp = gpar(fontsize = 8,fontface = "bold"),
                row_title = NULL,
                show_row_dend = TRUE,
                row_dend_width = unit(2, "cm"),
                column_dend_height = unit(2, "cm"),
                column_names_gp = gpar(fontsize = 9,fontface = "bold"),
                width = unit(8,"cm"))
write.table(file = paste0(data_add,"FullMatch_Heatmap_",length(filelist),"_",data_end,"_",reference,".txt"),
            mat,quote=FALSE,row.names=FALSE,sep="\t")
pdf(paste0(data_add,"FullMatch_Heatmap_",length(filelist),"_",data_end,"_",reference,".pdf"))
print(plt2)
dev.off()
