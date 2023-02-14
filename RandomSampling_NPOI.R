args=commandArgs(T) 
file1 = args[1] #  "nonpi6_sptd_filter.fa"
file2 = args[2] #  "pi6_sptd_filter.fa"
file3 = args[3] #   "1"
data_add = args[4]
cell_line = args[5]

#file1 = "nonpi6_sptd_filter.fa"
#file2 = "pi6_sptd_filter.fa"
#file3 = "1"
#data_add = "/home/user/Zenglin/degradome/datastore_final/sptd/"

print(file1)
A <- readLines(paste0(file1))
id1 <- A[seq(1,length(A),2)]
id2 <- A[seq(2,length(A),2)]
dfA <- data.frame("id"=id1, "sequence"=id2)

print(file2)
B <- readLines(file2)
id1 <- B[seq(1,length(B),2)]
id2 <- B[seq(2,length(B),2)]
dfB <- data.frame("id"=id1, "sequence"=id2)

select_row <- sample(1:nrow(dfA),nrow(dfB), replace = TRUE)
dfA1 <- dfA[select_row,]

output_file <- paste0(data_add,cell_line,"_nonPiRNA_",file3,".fa")
print(output_file)
idcol <- paste0(dfA1$id,"\n")
seqcol <- paste0(dfA1$sequence,"\n")
output <- paste0(idcol, seqcol)
writeLines(output,output_file,sep="")

