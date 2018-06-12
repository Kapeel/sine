library(plyr)
library(stringr)
library(rtracklayer)

#source('../CONFIG.R')

GENOME=commandArgs(trailingOnly=TRUE)
w=read.table(paste0(GENOME, '-matches.noTSD.MCSnames.8080.out'), header=F, sep='\t')


allw=read.table(paste0(GENOME, '-matches.noTSD.fa.fai'), header=F)
nofam=allw$V1[!allw$V1 %in% w$V1]
write(as.character(nofam), paste0(GENOME, '.RST.noExistingFam.txt'), ncolumns=1)
