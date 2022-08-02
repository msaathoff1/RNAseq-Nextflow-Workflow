#!/usr/bin/Rscript
# DESeq2 Script R-bit

args = commandArgs(trailingOnly = TRUE)

directory <- args[2]
format_table <- args[1]
sampleTable <- read.table(format_table, header = TRUE, sep = '\t')
sampleTable$condition <- factor(sampleTable$condition)

condition_a <- args[3]
condition_b <- args[4]

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)
dds <- DESeq(ddsHTSeq)
res <- results(dds)
write.csv(as.data.frame(res), file = paste("deseq2results/deseq2results_",condition_a,"_",condition_b,".csv",sep=""))

