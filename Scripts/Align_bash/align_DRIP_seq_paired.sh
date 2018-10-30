#!/bin/bash

#paths
quixpath="~/.guix-profile/bin"
bedGraphToBigWigpath="/home/kwreczy/programs"
bedtoolspath="/home/kwreczy/programs/bedtools-2.17.0/bin"
bbdukpath="/home/kwreczy/programs/bbmap"
fastqcpath="/home/kwreczy/programs/FastQC"
#bowtiepath="/home/kwreczy/programs/bowtie-1.1.2/"
akalindir="/data/akalin/Projects/AAkalin_HOTRegions/Data/DRIP_seq/hg19/"
hg19index="/data/akalin/Base/GenomeIndices/hg19/Bowtie2/hg19"
hg19chromInfo="/data/akalin/kasia/projects/HOT_DATA/hg19.chromInfo.txt"
adapter="/home/kwreczy/programs/bbmap/resources/truseq.fa.gz"

id=$1 #run id
dirname=$2

##mapping using Bowtie                                                                                   
# SINGLE-read                                                                                            
m=1 # number of allowed hits                                                                             
k=1 # number to report                                                                                   
n=1 #number of mismatches                                                             
procnum=3 # number of processors used       

#Log file from STDOUT for mapping
LOG2=$akalindir/$dirname/MAPPED/Bowtie/$id.log
#Error messages from STDERR for mapping
ERR2=$akalindir/$dirname/MAPPED/Bowtie/$id.err

~/.guix-profile/bin/bowtie2 --un-conc $akalindir/$dirname/MAPPED/Bowtie/output_file_un_conc --al-conc $akalindir/$dirname/MAPPED/Bowtie/output_file_al_conc -p 5 -x $hg19index -1 $akalindir/$dirname/RAW/fastq/$id'_1'.fastq -2 $akalindir/$dirname/RAW/fastq/$id'_2'.fastq  -S $akalindir/$dirname/MAPPED/Bowtie/$id.sam > $LOG2 > $ERR2







