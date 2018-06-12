#!/bin/bash -login
#set -e

#source ../CONFIG.sh

GENOME=$1
GENOMEFA=$2
SHORTID=$3
NEWSPECIES=$4

### the base path to find all the lcv alignments
HSDIR=/HelitronScanner
## where to find HelitronScanner.jar
HSJAR=/HelitronScanner/HelitronScanner/HelitronScanner.jar

## path to vsearch and silix for clustering
SILIX=silix
VSEARCH=/usr/local/vsearch-2.6.2-linux-x86_64/bin/vsearch

### helitron scanner needs some memory to load each chromosome in, so remember that when picking a queue
CPU=4
MEMGB=8

### Python version
PYTHON2=python2.7

### exsisting SINES
EXISTINGSINES=/vol_b/test_run2/sine/W22.RST.fa

echo Genome base name $GENOME
echo Genome fasta with path - needs to be up one from current dir: $GENOMEFA
echo CPUs: $CPU
echo Paths to silix and vsearch:  $SILIX, $VSEARCH
echo Path to python 2 executable: $PYTHON2

if [ $NEWSPECIES == 0 ]
then
echo Existing SINEs with families found : $EXISTINGSINES     ## need to actually check if file is present!
else
echo Generating all new family names for this species
fi

##################################
## Get Sine-Finder if necessary ##
##################################

if [ ! -f sine_finder.py ] 
then
#### get the SINE-Finder program from Wenke et al. 2011
wget http://www.plantcell.org/content/suppl/2011/08/29/tpc.111.088682.DC1/Supplemental_Data_Set_1-sine_finder.txt
#### change name
mv Supplemental_Data_Set_1-sine_finder.txt sine_finder.py
fi

##########################
#### run sinefinder   ####
##########################

echo Running SINE-Finder on $GENOME
#### I haven't been able to get sine_finder to work with reverse sequences, as it seems to report TSDs wrong on the reverse strand.
####   so I'm only reporting on the forward strand.
### -f both : outputs csv and fasta
$PYTHON2 sine_finder.py -T chunkwise -V1 -f both ${GENOMEFA}

#### sine_finder outputs the fasta with the TSD included. I remove these here, so they aren't considered when clustering into families
#mv ../${GENOME}-matches.fasta .
#mv ../${GENOME}-matches.csv .

$PYTHON2 remove_tsd_sinefinder.py ${GENOME}-matches.fasta ${GENOME}-matches.noTSD.fa

echo Clustering families
#################################################################################
if [ $NEWSPECIES == 1 ] ## We don't have existing 8 digit family names to add to, so we'll cluster everything
then
#### vsearch to identify homology, silix to cluster
$VSEARCH -allpairs_global ${GENOME}-matches.noTSD.fa -blast6out ${GENOME}-matches.noTSD.allvall.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads $CPU
# single linkage cluster those that are 80% identical to each other.
$SILIX ${GENOME}-matches.noTSD.fa ${GENOME}-matches.noTSD.allvall.8080.out -f SINE -i 0.8 -r 0.8 > ${GENOME}-matches.noTSD.8080.fnodes
## vsearch each against the MTEC exemplars

if [ ! -f TE_12-Feb-2015_15-35.fa ] ## get mtec db if needed
then
wget http://maizetedb.org/~maize/TE_12-Feb-2015_15-35.fa
fi
$VSEARCH --usearch_global TE_12-Feb-2015_15-35.fa -db ${GENOME}-matches.noTSD.fa -id 0.8 -query_cov 0.8 -target_cov 0.8 -blast6out ${GENOME}-matches.noTSD.TEDB8080.out -strand both -top_hits_only --threads $CPU

fi ## end loop if new species


#################################################################################
#### EXISTING FAMILIES   ########################################################
#################################################################################
if [ ! $NEWSPECIES == 1 ] ## we need to add to existing family names
then
#### vsearch to identify homology, silix to cluster - matching to EXISTING SINES!
$VSEARCH --usearch_global ${GENOME}-matches.noTSD.fa -db $EXISTINGSINES -id 0.8 -blast6out ${GENOME}-matches.noTSD.MCSnames.8080.out -query_cov 0.8 -target_cov 0.8 --threads $CPU --minseqlength 1

## index fasta so that we can subset entries by name
samtools faidx ${GENOME}-matches.noTSD.fa
#Rscript get_nofams.R ## this has hardcoded in to the r script w22 specific things!! 
Rscript get_nofams.R $GENOME
### get a fasta with each SINE that needs a new family assigned
xargs samtools faidx ${GENOME}-matches.noTSD.fa < ${GENOME}.RST.noExistingFam.txt >> ${GENOME}.RST.noExistingFam.fa

## generate new families for these, by clustering all-vs-all
${VSEARCH} -allpairs_global ${GENOME}.RST.noExistingFam.fa -blast6out ${GENOME}.RST.noExistingFam.allvall.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads 1
# single linkage cluster those that are 80% identical to each other.
$SILIX ${GENOME}.RST.noExistingFam.fa ${GENOME}.RST.noExistingFam.allvall.out -f SINE -i 0.8 -r 0.8 > ${GENOME}.RST.noExistingfam.8080.fnodes

## also get MTEC assignment for these new families
if [ ! -f TE_12-Feb-2015_15-35.fa ] ## get mtec db if needed
then
wget http://maizetedb.org/~maize/TE_12-Feb-2015_15-35.fa
fi
$VSEARCH --usearch_global TE_12-Feb-2015_15-35.fa -db ${GENOME}.RST.noExistingFam.fa -id 0.8 -query_cov 0.8 -target_cov 0.8 -blast6out ${GENOME}.RST.noExistingfam.TEDB8080.out -strand both -top_hits_only --threads $CPU

fi ## done with not new species


#######################################################################
### cluster into families and output final gff with this R script #####
#######################################################################

Rscript generate_gff_SINE.R $GENOME $GENOMEFA $SHORTID $NEWSPECIES

## and generate a fasta ready to go for the next annotator
$PYTHON2 switch_fasta_names.py ${GENOME}-matches.noTSD.fa ${GENOME}.RST.tabout > ${GENOME}.RST.fa
cat $EXISTINGSINES ${GENOME}.RST.fa > post-${GENOME}.existingRST.fa
