library(rtracklayer)
library(stringr)
library(data.table)
library(plyr)

#source('../CONFIG.R')

args=commandArgs(trailingOnly=TRUE)
GENOME=args[1]
GENOMEFA=args[2]
SHORTID=args[3]
NEWSPECIES=as.numeric(args[4])
#NEWSPECIES=args[4]

EXISTINGSINES <- c("/vol_b/test_run2/sine/W22.RST.fa")

print(paste0('Collating SINEs on ', GENOME, ' using short ID ', SHORTID))


if( NEWSPECIES == 0 ){
  print(paste0('Families exist in ', EXISTINGSINES, ', adding this genome to existing families'))
  
  ######### Read in silix results (family assignments to existing families)
  ## these match to existing families
  a=read.table(paste0(GENOME, '-matches.noTSD.MCSnames.8080.out'), header=F, stringsAsFactors=F)
  ## these are new families in this genome
  n=read.table(paste0(GENOME, '.RST.noExistingfam.8080.fnodes'), header=F, stringsAsFactors=F)
  ## find the maximum family number in existing families, increment.
  allSine=fread(paste0('grep ">" ', EXISTINGSINES), header=F)
  allSine$famNum=substr(allSine$V1,5,9)
  maxFamNum=as.numeric(max(allSine$famNum))
  
  
  ######### Read in sine_finder results, fix naming, classes
  sine=read.table(paste(GENOME, '-matches.csv', sep=''), header=T, stringsAsFactors=F)
  sine=unique(sine)   ## concerning - does it append to file? 
#  sine=sine[-which(sine$name=='name'),]    ## because this looks like a new header line!!
  sine$sequence=NULL
  sine$direct='+'
  sine$namecomp=paste(sine$name, sine$start, sine$end, 'F', paste('TSDlen', sine$TSD.len, sep=''), paste('TSDscore', sine$TSD.score, sep=''), paste('TSDmism', sine$TSD.mism, sep=''), sep='_')
  
  sine$existingFam=mapvalues(sine$namecomp, from=a$V1, to=substr(a$V2,1,8))
  sine$existingFam[grepl('TSD', sine$existingFam)]=NA
  
  sine$newFam=as.character(mapvalues(sine$namecomp, from=n$V2, to=as.character(n$V1)))
  sine$newFam[grepl('TSD', sine$newFam)]=NA

  ########### Add names to the new families
  ## first for new w22 families  ## this is super slow!!
  for (x in 1:length(table(sine$newFam))){
    famNum=str_pad(maxFamNum + x , 5, pad='0') # we want to increment families
    famName=names(rev(sort(table(sine$newFam))))[x]  ## and keep track of original families
    sine$newFamName[sine$newFam==famName & !is.na(sine$newFam)]=paste('RST', famNum, SHORTID, str_pad(1:sum(sine$newFam[!is.na(sine$newFam)]==famName), 5, pad='0'), sep='')
  }
  ### match new families to entries in the TE database (maizetedb.org)
  if(length(count.fields(paste0(GENOME, '.RST.noExistingfam.TEDB8080.out')))>0){ ## can't read in a file that is empty!
    pm=read.table(paste0(GENOME, '.RST.noExistingfam.TEDB8080.out'), header=F)
    pm$wicker3=str_split_fixed(pm$V1, '_', 3)[,1]
    pm=pm[pm$wicker3=='RST',]
    pm$family=str_split_fixed(pm$V1, '_', 3)[,2]
    
    pm.gd=merge(pm, a, by.x='V2', by.y='V2', all.x=T)
  
    sine$mtec.fam=NA
    for (i in 1:nrow(pm.gd)){
            sine$mtec.fam[sine$V1==pm.gd$V1.y[i]]=pm.gd$family[i]
            }
  } else{ sine$mtec.fam=NA }
  
  ## assign 11 digit copy name (SHORTIDXXXXX) for each copy in an existing family
  sine$Name=NA
  for (x in names(table(sine$existingFam))){
    sine$Name[sine$existingFam==x & !is.na(sine$existingFam)]=paste(x, SHORTID, str_pad(1:sum(sine$existingFam[!is.na(sine$existingFam)]==x), 5, pad='0'), sep='')
  }
  sine$Name[is.na(sine$existingFam)]=sine$newFamName[is.na(sine$existingFam)]
  
  ## done!!!!
  print(paste0('There are ', sum(is.na(sine$Name)), ' SINEs without names'))
} ## end clustering with existing families

if( NEWSPECIES == 1){
  print(paste0('Making new families for', GENOME))
  ######### Read in silix results (family assignments)
  a=read.table(paste(GENOME, '-matches.noTSD.8080.fnodes', sep=''))
  
  ### Switch SINE to RST, the wicker 3 letter code
  ## because these are plants, and found with the rna pol A and B boxes, these are all RST
  ## remove one digit because there are not 100,000 families, one order of mag less
  #f$V1=sub('SINE0', 'RST', f$V1)  ### let's switch and do as early as we can!!!
  a$V1=sub('SINE0', 'RST', a$V1)
  ######### Read in sine_finder results, fix naming, classes
  sine=read.table(paste(GENOME, '-matches.csv', sep=''), header=T)
  sine$sequence=NULL
  sine$direct='+'
  sine$namecomp=paste(sine$name, sine$start, sine$end, 'F', paste('TSDlen', sine$TSD.len, sep=''), paste('TSDscore', sine$TSD.score, sep=''), paste('TSDmism', sine$TSD.mism, sep=''), sep='_')
  
  f=merge(sine, a, by.x='namecomp', by.y='V2', all=T)
  f$fam=table(f$V1)[f$V1] 
  ### give each TE copy a unique identifier
  f$Name=''
  ### can't figure out how to assign properly, so giving up and doing a for loop
  #f$Name=sapply(names(table(f$V1)), function(x) f$Name[f$V1==x]=paste(x, 'B73v4', str_pad(1:sum(f$V1==x), 5, pad='0'), sep=''))
  for (x in names(table(f$V1))){
          f$Name[f$V1==x]=paste(x, SHORTID, str_pad(1:sum(f$V1==x), 5, pad='0'), sep='')
          }
          
  
  
  ### match to entries in the TE database (maizetedb.org)
  pm=read.table(paste(GENOME, '-matches.noTSD.TEDB8080.out', sep=''), header=F)
  pm$wicker3=str_split_fixed(pm$V1, '_', 3)[,1]
  pm=pm[pm$wicker3=='RST',]
  #pm.gd=f$namecomp %in% pm$V2
  pm$family=str_split_fixed(pm$V1, '_', 3)[,2]
  
  pm.gd=merge(pm, a, by.x='V2', by.y='V2', all.x=T)

  f$mtec.fam=NA
  for (i in 1:nrow(pm.gd)){
          f$mtec.fam[f$V1==pm.gd$V1.y[i]]=pm.gd$family[i]
          }
  sine=f
} ## close new species



####################
## output the gff ##
####################

sine.gff=data.frame(sine$name, 'SineFinder', 'SINE_element', sine$start, sine$end, '.', '+', '.', paste('ID=', sine$Name, sep=''))
write.table(sine.gff, paste(GENOME, '.RST.gff3', sep=''), quote=F, sep='\t', row.names=F, col.names=F)

### also keep track of fasta names and newly assigned gffnames as to easily convert between the two (e.g. switching fasta names)
sine.out=data.frame(TEID=sine$Name, namecompare=sine$namecomp, TSD.len=sine$TSD.len, TSD.score=sine$TSD.score, TSD.mism=sine$TSD.mism, mtec.fam=sine$mtec.fam)
write.table(sine.out, paste(GENOME, '.RST.tabout', sep=''), quote=F, sep='\t', row.names=F, col.names=T)





