#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=08:0:00

# 20.02.2017

# paths
quixpath="~/.guix-profile/bin"
bedGraphToBigWigpath="/home/kwreczy/programs"
bedtoolspath="/home/kwreczy/programs/bedtools-2.17.0/bin"
bbdukpath="/home/kwreczy/programs/bbmap"
fastqcpath="/home/kwreczy/programs/FastQC"

chromInfo="/data/akalin/Base/GenomeIndices/ce10/WholeGenomeFasta/chrNameLength.txt" 
INDEX="/data/akalin/Base/GenomeIndices/ce10/Bowtie2/ce10"

# number of threads
NMBR_THREADS=20

# http://bioinformatics.mdc-berlin.de/intro2UnixandSGE/sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html
RUN_QSUB="qsub -V -l h_vmem=5G -pe smp 15 "

# samples
DIR="/data/akalin/Projects/AAkalin_HOTRegions/Data/DRIP_seq/ce10/"
cd $DIR

RUNS=(SRR3993994 SRR3993996 SRR3993998)


# FASTQC for paired-end
for id in ${RUNS[*]};
do
$fastqcpath/fastqc --threads $NMBR_THREADS -o $DIR/RAW/QC  $DIR/RAW/fastq/$id'_1'.fastq  2> $DIR/RAW/QC/$id'_1'.log
$fastqcpath/fastqc --threads $NMBR_THREADS -o $DIR/RAW/QC  $DIR/RAW/fastq/$id'_2'.fastq  2> $DIR/RAW/QC/$id'_2'.log
done

# bbduk2 for paired-end
for id in ${RUNS[*]};
do
 /home/kwreczy/programs/bbmap/bbduk2.sh -Xmx1g  in1=$DIR/RAW/fastq/$id'_1'.fastq in2=$DIR/RAW/fastq/$id'_2'.fastq out1=$DIR/Cleaned/$id'_1'.bbduk.fastq out2=$DIR/Cleaned/$id'_2'.bbduk.fastq t=10 qtrim=r trimq=36 hdist=1 minlength=18 overwrite=true k=25  2> $DIR/Cleaned/$id.bbduk.log
done

# FASTQC after bbduk2 for paired-end
for id in ${RUNS[*]};
do
$fastqcpath/fastqc --threads $NMBR_THREADS -o $DIR/Cleaned/QC $DIR/Cleaned/$id'_1'.bbduk.fastq 2> $DIR/Cleaned/QC/$id'_1'.bbduk.log
$fastqcpath/fastqc --threads $NMBR_THREADS -o $DIR/Cleaned/QC  $DIR/Cleaned/$id'_2'.bbduk.fastq  2> $DIR/Cleaned/QC/$id'_2'.bbduk.log
done

# bowtie2 for paired-end
for id in ${RUNS[*]};
do
echo $id
LOG_bowtie2=$DIR/Mapped/$id.bbduk.bowtie2.log
un_conc=$DIR/Mapped/$id.bowtie2.output_file_un_conc
al_conc=$DIR/Mapped/$id.bowtie2.output_file_al_conc
in1=$DIR/Cleaned/$id'_1'.bbduk.fastq
in2=$DIR/Cleaned/$id'_2'.bbduk.fastq
out=$DIR/Mapped/$id.bbduk.bowtie2.sam
 /home/kwreczy/.guix-profile/bin/bowtie2 -x $INDEX -1 $in1 -2 $in2 -S $out --un-conc $un_conc --al-conc $al_conc -p $NMBR_THREADS 2> $LOG_bowtie2

done

# convert SAM to BAM, sort and index BAM
for id in ${RUNS[*]};
do
samtools view -Sb  $DIR/Mapped/$id.bbduk.bowtie2.sam  >  $DIR/Mapped/$id.bbduk.bowtie2.bam
samtools sort $DIR/Mapped/$id.bbduk.bowtie2.bam  $DIR/Mapped/$id.bbduk.bowtie2.sorted
samtools index $DIR/Mapped/$id.bbduk.bowtie2.sorted.bam
done


# convert BAM to BigWig
for id in ${RUNS[*]}; do

# BAM -> BED
$bedtoolspath/bamToBed -i $DIR/Mapped/$id.bbduk.bowtie2.sorted.bam > $DIR/Mapped/$id.bed
# BED -> bedGRaph
$bedtoolspath/bedtools genomecov -bg -i $DIR/Mapped/$id.bed -g $chromInfo > $DIR/Mapped/$id.bedgraph
sort -k1,1 -k2,2n $DIR/Mapped/$id.bedgraph > $DIR/Mapped/$id.sorted.bedgraph
# BED -> BigWig
$bedGraphToBigWigpath/bedGraphToBigWig $DIR/Mapped/$id.sorted.bedgraph $chromInfo $DIR/Mapped/$id.bw

done



