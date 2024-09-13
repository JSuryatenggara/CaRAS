# Load required libraries
pacs <- c("DESeq2")

if(sum(as.numeric(!pacs %in% installed.packages())) != 0){
  ins <- pacs[!pacs %in% installed.packages()]
  for(i in 1:length(ins)) {
    install.packages(ins, dependencies = T)
    break()}
  sapply(pacs, require, character = T) 
} else {
  sapply(pacs, require, character = T) 
}


library(DESeq2)

args = commandArgs(trailingOnly=TRUE)

input = args[1]
output = args[2]

counts = read.delim(file = input,header=T,row.names = 1,check.names = F)
sf = estimateSizeFactorsForMatrix(counts)
sf = as.data.frame(sf)
sf$normalizer = (1/sf$sf)

write.table(sf,file=output,quote=F,sep='\t',col.names=T)
