#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=08:0:00

# 24.02.2017

#paths
quixpath="~/.guix-profile/bin"
bedGraphToBigWigpath="/home/kwreczy/programs"
bedtoolspath="/home/kwreczy/programs/bedtools-2.17.0/bin"
bbdukpath="/home/kwreczy/programs/bbmap"
fastqcpath="/home/kwreczy/programs/FastQC"
INDEX="/data/akalin/Base/GenomeIndices/mm9/Bowtie2/mm9"
chromInfo="/home/kwreczy/mm9.chromInfo.txt"
adapter="/home/kwreczy/programs/bbmap/resources/adapters.fa" #Illumina Universal adapter are there

# number of threads
NMBR_THREADS=20

# command for cluster
RUN_QSUB="qsub -V -l h_vmem=5G -pe smp 15 "

# samples
DIR="/data/akalin/Projects/AAkalin_HOTRegions/Data/DRIP_seq/mm9/"
cd $DIR
RUNS=(SRR1952485 SRR1952486 SRR1952450)

## TODO add line to gunzip files if there are compressed


# echo "FASTQC"
# # FASTQC for single-end
# mkdir $DIR/RAW/QC
# for id in ${RUNS[*]};
# do
# echo $id
# $fastqcpath/fastqc --threads $NMBR_THREADS -o $DIR/RAW/QC  $DIR/RAW/fastq/$id.fastq  2> $DIR/RAW/QC/$id.log
# done

mkdir $DIR/Cleaned
for id in ${RUNS[*]};
do
echo $id
 /home/kwreczy/programs/bbmap/bbduk2.sh -Xmx1g  in=$DIR/RAW/fastq/$id.fastq out=$DIR/Cleaned/$id.bbduk.fastq  t=10 qtrim=r trimq=28 hdist=1 minlength=18 overwrite=true k=25  2> $DIR/Cleaned/$id.bbduk.log
done

echo "FASTQC after bbduk2 for single-end"
# FASTQC after bbduk2 for single-end
mkdir $DIR/Cleaned/QC
for id in ${RUNS[*]};
do
echo $id
$fastqcpath/fastqc --threads $NMBR_THREADS -o $DIR/Cleaned/QC $DIR/Cleaned/$id.bbduk.fastq 2> $DIR/Cleaned/QC/$id.bbduk.log
done

echo "single-end"
# bowtie2 for single-end
mkdir $DIR/Mapped
for id in ${RUNS[*]};
do
echo $id
LOG_bowtie2=$DIR/Mapped/$id.bbduk.bowtie2.log
un_conc=$DIR/Mapped/$id.bowtie2.output_file_un_conc
al_conc=$DIR/Mapped/$id.bowtie2.output_file_al_conc
in=$DIR/Cleaned/$id.bbduk.fastq
out=$DIR/Mapped/$id.bbduk.bowtie2.sam
/home/kwreczy/.guix-profile/bin/bowtie2 -x $INDEX -U $in -S $out --un-conc $un_conc --al-conc $al_conc -p $NMBR_THREADS 2> $LOG_bowtie2

done

echo "convert SAM to BAM, sort and index BAM"
# convert SAM to BAM, sort and index BAM
for id in ${RUNS[*]};
do
/home/kwreczy/.guix-profile/bin/samtools view -Sb  $DIR/Mapped/$id.bbduk.bowtie2.sam  >  $DIR/Mapped/$id.bbduk.bowtie2.bam
/home/kwreczy/.guix-profile/bin/samtools sort $DIR/Mapped/$id.bbduk.bowtie2.bam  $DIR/Mapped/$id.bbduk.bowtie2.sorted
/home/kwreczy/.guix-profile/bin/samtools index $DIR/Mapped/$id.bbduk.bowtie2.sorted.bam
done


