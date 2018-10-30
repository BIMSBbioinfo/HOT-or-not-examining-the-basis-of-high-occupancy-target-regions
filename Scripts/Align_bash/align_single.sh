#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=9G




#paths
quixpath="~/.guix-profile/bin"
bedGraphToBigWigpath="/home/kwreczy/programs"
bedtoolspath="/home/kwreczy/programs/bedtools-2.17.0/bin"
bbdukpath="/home/kwreczy/programs/bbmap"
fastqcpath="/home/kwreczy/programs/FastQC"
bowtiepath="/home/kwreczy/programs/bowtie-1.1.2/"


######
akalindir="/data/akalin/Projects/AAkalin_HOTRegions/Data/DRIP_seq/hg19/"


hg19index="/data/akalin/Base/GenomeIndices/hg19/Bowtie/hg19"
hg19chromInfo="/data/akalin/kasia/projects/HOT_DATA/hg19.chromInfo.txt"
adapter="/home/kwreczy/programs/bbmap/resources/truseq.fa.gz"

id=$1 #run id
dirname=$2

echo $id
echo $dirname
echo $akalindir/$dirname/RAW/fastq/$id.fastq.gz

##QC on fastq files                                                                                      
echo "QC"
                         
$fastqcpath/fastqc -o $akalindir/$dirname/RAW/QC  $akalindir/$dirname/RAW/fastq/$id.fastq.gz  2> $akalindir/$dirname/RAW/QC/$id.log                                                                                                                               


##bbduk2 - Adapter/Quality Trimming and Filtering                                                                                  
# ktrim=f  Trim reads to remove bases matching reference kmers. r (trim to the right),                                             
# kmask=f  Replace bases matching ref kmers with another symbol.                                                                   
#          Allows any non-whitespace character other than t or f,                                                                  
#          and processes short kmers on both ends.  'kmask=lc' will                                                                
#          convert masked bases to lowercase.                                                                                      
# mink=0   Look for shorter kmers at read tips down to this length,                                                                
#         when k-trimming or masking.  0 means disabled.  Enabling                                                                 
#         this will disable maskmiddle.                                                                                            
# qtrim=f     Trim read ends to remove bases with quality below trimq.Performed AFTER looking for kmers.rl (trim both ends),       
# trimq=6     Regions with average quality BELOW this will be trimmed.                                                             
# minlength   Reads shorter than this after trimming will be discarded.  Pairs will be discarded if both are shorter.              
# k=27         Kmer length used for finding contaminants.  Contaminants                                                            
# shorter than k will not be found.  k must be at least 1.                                                                         
$bbdukpath/bbduk2.sh -Xmx1g qin=33 in=$akalindir/$dirname/RAW/fastq/$id.fastq.gz out=$akalindir/$dirname/RAW/cleaned/bbduk2/$id.qin33.fastq minlength=18 qtrim=r trimq=20 ktrim=r k=25 mink=11 ref=$adapter hdist=1 overwrite=true   2> $akalindir/$dirname/RAW/cleaned/bbduk2/$id.qin=33.log
$bbdukpath/bbduk2.sh -Xmx1g ignorequality in=$akalindir/$dirname/RAW/fastq/$id.fastq.gz out=$akalindir/$dirname/RAW/cleaned/bbduk2/$id.ignorequality.fastq minlength=18 qtrim=r trimq=20 ktrim=r k=25 mink=11 ref=$adapter hdist=1 overwrite=true   2> $akalindir/$dirname/RAW/cleaned/bbduk2/$id.ignorequality.log

cp $akalindir/$dirname/RAW/cleaned/bbduk2/$id.qin33.fastq $akalindir/$dirname/RAW/cleaned/bbduk2/$id.fastq


#fastqc after bbduk2
$fastqcpath/fastqc -o $akalindir/$dirname/RAW/cleaned/QC  $akalindir/$dirname/RAW/cleaned/bbduk2/$id.fastq  2> $akalindir/$dirname/RAW/cleaned/QC/$id.log


echo "bowtie"

##mapping using Bowtie                                                                                   
# SINGLE-read                                                                                            
m=1 # number of allowed hits                                                                             
k=1 # number to report                                                                                   
n=1 #number of mismatches                                                             
procnum=3 # number of processors used                                                                                                                                                            

$bowtiepath/bowtie -m $m -S -k $k -n $n -p $procnum  $hg19index $akalindir/$dirname/RAW/cleaned/bbduk2/$id.fastq  > $akalindir/$dirname/MAPPED/Bowtie/$id.sam 2> $akalindir/$dirname/MAPPED/Bowtie/$id.log


##SAM to BAM format                                                                                                                                                                              

echo "sam2bam"

samtools view -Sb  $akalindir/$dirname/MAPPED/Bowtie/$id.sam  >  $akalindir/$dirname/MAPPED/Bowtie/$id.bam


##sort and index alignment                                                                                                  
samtools sort $akalindir/$dirname/MAPPED/Bowtie/$id.bam  $akalindir/$dirname/MAPPED/Bowtie/$id.sorted

samtools index $akalindir/$dirname/MAPPED/Bowtie/$id.sorted.bam



##BAM -> BED                                                                                                                                        
$bedtoolspath/bamToBed -i $akalindir/$dirname/MAPPED/Bowtie/$id.sorted.bam > $akalindir/$dirname/MAPPED/$id.bed

##BED -> bedGRaph                                                                                                                                   
$bedtoolspath/bedtools genomecov -bg -i $akalindir/$dirname/MAPPED/$id.bed -g $hg19chromInfo > $akalindir/$dirname/MAPPED/$id.bedgraph
sort -k1,1 -k2,2n $akalindir/$dirname/MAPPED/$id.bedgraph > $akalindir/$dirname/MAPPED/$id.sorted.bedgraph

echo "bed2bigwig"

##BED -> BigWig                                                                                                                                     
$bedGraphToBigWigpath/bedGraphToBigWig $akalindir/$dirname/MAPPED/$id.sorted.bedgraph $hg19chromInfo $akalindir/$dirname/MAPPED/BigWig/$id.bw




