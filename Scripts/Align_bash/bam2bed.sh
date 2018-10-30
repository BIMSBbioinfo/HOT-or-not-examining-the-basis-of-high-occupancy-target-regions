#!/bin/bash
#$ -pe smp 5
#$ -l h_vmem=10G

#paths
quixpath="~/.guix-profile/bin"
bedGraphToBigWigpath="/home/kwreczy/programs"
bedtoolspath="/home/kwreczy/programs/bedtools-2.17.0/bin"

id=$1 # sample id
dirname=$2 # path to directory with sam file
chromInfo=$3 # e.g. path to infor about chrom length hg19chromInfo="/data/akalin/kwreczy/projects/HOT_DATA/hg19.chromInfo.txt"

# Convert sam to bam
samtools view -Sb  $dirname/$id.sam  >  $dirname/$id.bam

# Sort bam file
# It adds ".bam" suffix to file $id.sorted automatically
samtools sort $dirname/$id.bam  $dirname/$id.sorted

# Index bam file
samtools index $dirname/$id.sorted.bam

##BAM -> BED                                                                                                                                        
$bedtoolspath/bamToBed -i $dirname/$id.sorted.bam > $dirname/$id.bed

##BED -> bedGRaph                                                                                                                                   
$bedtoolspath/bedtools genomecov -bg -i $dirname/$id.bed -g $chromInfo > $dirname/$id.bedgraph
sort -k1,1 -k2,2n $dirname/$id.bedgraph > $dirname/$id.sorted.bedgraph

##BED -> BigWig                                                                                                                                     
$bedGraphToBigWigpath/bedGraphToBigWig $dirname/$id.sorted.bedgraph $chromInfo $dirname/$id.bw




