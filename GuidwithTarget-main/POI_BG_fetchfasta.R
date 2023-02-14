args=commandArgs(T) 
file1=args[1]
piRNA=args[2]
piRNA_prefix <- unlist(strsplit(file1,"_"))[1]
#file1 <- "ago3_filter_cnt.sam"
#piRNA <- "../pidroso/piRNAdb.dme.v1_7_6.fa"

library("data.table")
input <- fread(file1)
input.stat <- as.data.frame(table(input$V1))

`%notin%` <- Negate(`%in%`)

Fasta2df2 <- function(filename){
             allfasta  = readLines(filename)
             id1 <- allfasta[seq(1,length(allfasta),2)]
             id2 <- allfasta[seq(2,length(allfasta),2)]
             df <- data.frame("id"=id1, "seq"=id2)
             return(df)
			 }

Df2fasta <- function(df2,filename){
		idcol <- paste0(df2$id,"\n")
		seqcol <- paste0(df2$seq,"\n")
		output <- paste0(idcol,seqcol)
		writeLines(output,filename,sep="")
}
allpi <- Fasta2df2(piRNA)
POI <- allpi[which(allpi$id %in% paste0(">",input.stat$Var1)),]
NPOI <- allpi[which(allpi$id %notin% paste0(">",input.stat$Var1)),]
Df2fasta(POI,paste0(piRNA_prefix,"_IP.fa"))
Df2fasta(NPOI,paste0(piRNA_prefix,"_BG.fa"))
print("Finished!")
