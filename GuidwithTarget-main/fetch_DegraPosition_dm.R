args=commandArgs(T)
file1 = args[1] # wtago3_candidate.sam
file2 = args[2] # mutago3_candidate.sam
data_add = args[3] # address:/data3/zenglin/temp/LinZeng/degradome/NGdata/"
reference = args[4] # reference

options(scipen = 100)
library("data.table")
suppressPackageStartupMessages(library("dplyr"))


cellline = unlist(strsplit(file1,"_"))[1]
cellline = substr(cellline,3,nchar(cellline))

flagstat1 <- readLines(paste0(data_add,"wt",cellline,"_flagstat.txt"))
flagstat_wt <- as.numeric(unlist(strsplit(flagstat1[5]," "))[1])

flagstat2 <- readLines(paste0(data_add,"mut",cellline,"_flagstat.txt"))
flagstat_mut <- as.numeric(unlist(strsplit(flagstat2[5]," "))[1])

########################Part1
print("Processing wide type:")
wt <- fread(paste0(data_add,file1),sep="\t")
col_2_rename = c("id","strand","chr","pos_1_based")
colnames(wt) <- col_2_rename
wt$strand <- ifelse(wt$strand==0,"+","-")  #!!!!!!!!!!!!!!!!!!!!!!!!!!!bed file 0-based  not include last character
range2ext = 50
if(reference=="1"){
  wt <- wt[which(wt$strand=="+"),]
  range2ext = 30
}
wt$start <- as.numeric(wt$pos_1_based) - range2ext - 1
wt$start <- ifelse(wt$start<0,0,wt$start)
wt$end <- as.numeric(wt$pos_1_based) + range2ext
wt$id <- unlist(strsplit(wt$id,"\\."))[seq(2,nrow(wt)*2,2)]
wt$uniqpos <- paste0(wt$chr,"_",wt$pos_1_based,"_",wt$strand)
products_stat <- as.data.frame(table(wt$uniqpos)) 
wt <- merge(wt,products_stat,by.x="uniqpos",by.y="Var1")
wt <- wt %>% select(uniqpos,id,strand,chr,pos_1_based,start,end,Freq) %>% group_by(uniqpos) %>% mutate(combi_id = paste(id, collapse = " | "))
wt <- unique(as.data.frame(wt[,c("chr","start","end","uniqpos","Freq","strand","combi_id","pos_1_based")]))
write.table(file=paste0(data_add,unlist(strsplit(file1,"_"))[1],"_3CleavProcutsStatRaw.txt"),wt,sep="\t",quote=FALSE,row.names=FALSE)
reduced_wt <- wt[,c("chr","start","end","uniqpos","Freq","strand","pos_1_based")]
reduced_wt$Freq <- round((as.numeric(reduced_wt$Freq) + 1)/flagstat_wt * 1000000,4)
reduced_wt <- reduced_wt[order(reduced_wt$Freq,decreasing=TRUE),]
colnames(reduced_wt) = c("chr","start","end","uniqpos","FreqWT","strand","pos_1_based")
write.table(file=paste0(data_add,unlist(strsplit(file1,"_"))[1],"_3CleavProductsStatNorm.txt"),reduced_wt,sep="\t",quote=FALSE,row.names=FALSE)
print("wide type processed!")
rm(wt)


#########################Part2
print("Processing Mut:")
mut <- fread(paste0(data_add,file2),sep="\t")
col_2_rename = c("id","strand","chr","pos_1_based")
colnames(mut) <- col_2_rename

mut$strand <- ifelse(mut$strand==0,"+","-")
if(reference=="1"){
  mut <- mut[which(mut$strand=="+"),]
}
mut$start <- as.numeric(mut$pos_1_based) - range2ext - 1
mut$start <- ifelse(mut$start<0,0,mut$start)
mut$end <- as.numeric(mut$pos_1_based) + range2ext
mut$id <- unlist(strsplit(mut$id,"\\."))[seq(2,nrow(mut)*2,2)]
mut$uniqpos <- paste0(mut$chr,"_",mut$pos_1_based,"_",mut$strand)
products_stat <- as.data.frame(table(mut$uniqpos))
mut <- merge(mut,products_stat,by.x="uniqpos",by.y="Var1")
mut <- mut %>% select(uniqpos,id,strand,chr,pos_1_based,start,end,Freq) %>% group_by(uniqpos) %>% mutate(combi_id = paste(id, collapse = " | "))
mut <- unique(as.data.frame(mut[,c("chr","start","end","uniqpos","Freq","strand","combi_id","pos_1_based")]))
write.table(file=paste0(data_add,unlist(strsplit(file2,"_"))[1],"_3CleavProcutsStatRaw.txt"),mut,sep="\t",quote=FALSE,row.names=FALSE)
reduced_mut <- mut[,c("chr","start","end","uniqpos","Freq","strand","pos_1_based")]
reduced_mut$Freq = round((as.numeric(reduced_mut$Freq) + 1)/flagstat_mut * 1000000,4) 
reduced_mut <- reduced_mut[order(reduced_mut$Freq,decreasing=TRUE),]
colnames(reduced_mut) = c("chr","start","end","uniqpos","FreqMut","strand","pos_1_based")
write.table(file=paste0(data_add,unlist(strsplit(file2,"_"))[1],"_3CleavProductsStatNorm.txt"),reduced_mut,sep="\t",quote=FALSE,row.names=FALSE)
rm(mut)

##########################Part3
reduced_mut <- merge(reduced_mut, reduced_wt[,c("uniqpos","FreqWT")],by="uniqpos",all.x=TRUE)
reduced_wt <- merge(reduced_wt,reduced_mut[,c("uniqpos","FreqMut")],by="uniqpos",all.x=TRUE)
reduced_total <- unique(rbind(reduced_mut,reduced_wt))

#temp_total[is.na(temp_total)] <- 0 #ERROR
reduced_total[which(is.na(reduced_total$FreqWT)),"FreqWT"] = round(1/flagstat_wt * 1000000,4)
reduced_total[which(is.na(reduced_total$FreqMut)),"FreqMut"] = round(1/flagstat_mut * 1000000,4)
reduced_total$FreqComp <- paste0(reduced_total$FreqMut,"|",reduced_total$FreqWT,"|",reduced_total$pos_1_based)
reduced_total$placeholder = 255

reduced_total$start <- as.numeric(reduced_total$start)
reduced_total$end <- as.numeric(reduced_total$end)
reduced_total <- reduced_total[which(reduced_total$FreqMut < reduced_total$FreqWT),]
reduced_total <- reduced_total[order(reduced_total$FreqWT,decreasing=TRUE),] #1345442

print("mut file processed!")

write.table(file=paste0(data_add,cellline,"_insert.txt"),
  reduced_total,
  sep="\t",
  quote=FALSE,
  row.names=FALSE)
write.table(file=paste0(data_add,cellline,"_insert.bed"),
  reduced_total[,c("chr","start","end","FreqComp","placeholder","strand")],
  sep="\t",
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE)





