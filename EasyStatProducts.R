args=commandArgs(T)
file1 = args[1]
#file = "ago3_Parsed_Filter_POI.txt"
result <- read.table(file1,sep="@")
names2sign <- c("query_id","sub_id","q_start","q_end","s_tart","s_end",
              "q_match_pos","q_mismatch_pos","alignment_length",
              "match_num","mismatch_num","identity","dis1",
              "pos1_based","start_abs","dis","relpos")
colnames(result) <- names2sign
extractpos <- function(matchpos){
  return (length(unlist(strsplit(matchpos,","))))
}

result$deg_wt = as.numeric(unlist(strsplit(result$sub_id,"\\|"))[seq(2,nrow(result)*3,3)])
result$mismatch_length <- unlist(lapply(result$q_mismatch_pos,extractpos))
print(paste0("0.Total Cleaved 3 Products decreased in Mut Compared WT: ",sum(result$deg_wt)))
print(paste0("  Total Records : ",nrow(result)," records!"))

fullmatch <- result[which(result$identity=="100%"),]
print(paste0("1.Explianed by fullmatch : ",sum(fullmatch$deg_wt),"(",
             round((sum(fullmatch$deg_wt)/sum(result$deg_wt)),4)*100,"%)"))
print(paste0("  Explianed by fullmatch : ",nrow(fullmatch)," records!"))

singlemiss <- result[which(result$mismatch_length=="1"),]
print(paste0("2.Explianed by single miss : ",
             sum(singlemiss$deg_wt),
             "(",
             round((sum(singlemiss$deg_wt)/sum(result$deg_wt)),4)*100,
             "%)"))


doublemiss <- result[which(result$mismatch_length=="2"),]
print(paste0("2.Explianed by double miss : ",
             sum(doublemiss$deg_wt),
             "(",
             round((sum(doublemiss$deg_wt)/sum(result$deg_wt)),4)*100,
             "%)"))
triplemiss <- result[which(result$mismatch_length=="3"),]
print(paste0("3.Explianed by triple miss : ",
             sum(triplemiss$deg_wt),
             "(",
             round((sum(triplemiss$deg_wt)/sum(result$deg_wt)),4)*100,
             "%)"))
print("==============================================================")


#############################################################
fullmatch <- result[which(result$identity=="100%"),]
print(paste0("Explianed by fullmatch : ",sum(fullmatch$deg_wt),"(",
             round((sum(fullmatch$deg_wt)/sum(result$deg_wt)),4)*100,"%)"))
for(start in sort(unique(fullmatch$q_start))){
  startfullmatch <- fullmatch[which(fullmatch$q_start==start),]
  print(paste0("      Explianed by fullmatch G",start,": ",
               sum(startfullmatch$deg_wt),
               "(",
               round((sum(startfullmatch$deg_wt)/sum(result$deg_wt)),4)*100,
               "%)"))
  for(end in sort(unique(startfullmatch$q_end))){
    endfullmatch <- startfullmatch[which(startfullmatch$q_end==end),]
    print(paste0("        Explianed by fullmatch G",start,"-","G",end," ",
                 sum(endfullmatch$deg_wt),
                 "(",
                 round((sum(endfullmatch$deg_wt)/sum(result$deg_wt)),4)*100,
                 "%)"))
  }
}

print("======================================================")

#############################################################
singlemiss <- result[which(result$mismatch_length=="1"),]
print(paste0("    Explianed by singlemiss : ",
             sum(singlemiss$deg_wt),
             "(",
             round((sum(singlemiss$deg_wt)/sum(result$deg_wt)),4)*100,
      "%)"))
for(start in sort(unique(singlemiss$q_start))){
  startmiss <- singlemiss[which(singlemiss$q_start==start),]
  print(paste0("     Explianed by single miss G",start,": ",
               sum(startmiss$deg_wt),
               "(",
               round((sum(startmiss$deg_wt)/sum(result$deg_wt)),4)*100,
               "%)"))
  startmiss$q_end <- as.numeric(startmiss$q_end)
  for(end in sort(unique(startmiss$q_end))){
    endmiss <- startmiss[which(startmiss$q_end==end),]
    endmiss$q_mismatch_pos <- as.numeric(endmiss$q_mismatch_pos)
      for(snp in sort(unique(endmiss$q_mismatch_pos))){
        snpmiss <- endmiss[which(endmiss$q_mismatch_pos==snp),]
        print(paste0("        Explianed by single miss G",start,"-","G",end,"-",
                     "M",snp," ",
                     sum(snpmiss$deg_wt),
                     "(",
                     round((sum(snpmiss$deg_wt)/sum(result$deg_wt)),4)*100,
                     "%)"))
      }
  }
}



