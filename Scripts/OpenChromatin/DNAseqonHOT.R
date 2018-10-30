
options(scipen=999)

library(ggplot2)
library(genomation)
library(Rsamtools)
library(ggplot2)
library("gridExtra")
library(reshape)
library(cowplot)
library(plyr)
library(reshape)


source("~/projects/HOT-regions/basic_functions.R")
library(BS.genomeSpecies, character.only = TRUE)
add.intervals = function(hot){
  hot$x = "0"
  hot[hot$perc >= 0.99 &  hot$perc <=1]$x  = "0.99-1"
  hot[hot$perc >= 0.9 &  hot$perc  <0.99]$x = "0.9-0.99"
  hot[hot$perc >= 0.8 &  hot$perc <0.9]$x = "0.8-0.9"
  hot[hot$perc >= 0.7 &  hot$perc <0.8]$x = "0.7-0.8"
  hot[hot$perc >= 0.6 &  hot$perc <0.7]$x = "0.6-0.7"
  hot[hot$perc >= 0.5 &  hot$perc <0.6]$x = "0.5-0.6"
  hot[hot$perc >= 0.4 &  hot$perc <0.5]$x = "0.4-0.5"
  hot[hot$perc >= 0  &  hot$perc < 0.4]$x = "0-0.4"
  hot$x = as.factor(hot$x)
  return(hot)
}

#############################################################################
## hg19
#############################################################################


hg19 = processHOTsummits_ucsc("/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/HOThg19_filterNearBy2kb_Uniform_nopol2a_summits.bed",
                              2000, "hg19")
# remove unnecessary chromosomes
hg19 = keepSeqlevels(hg19,
                     unique(as.character(seqnames(hg19))) )
black = readBed("/home/buyar/datasets/genomes/hg19/wgEncodeDacMapabilityConsensusExcludable.bed")
hg19noblack = hg19[which(hg19 %outside% black)]
hg19noblack = hg19noblack[order(-hg19noblack$perc)]
hot=hg19noblack
wind = hg19noblack
wind = add.intervals(wind)
wind = keepStandardChromosomes(wind, "Homo_sapiens", pruning.mode="tidy")
wind.inter = split(wind, wind$x)

#HOT (>99th percentile), MILD (between 99th and 75th percentile), and COLD regions (below 75th percentile)
HOT=wind[ wind$perc>.99 ]
MILD=wind[ wind$perc<.99 & wind$perc>.75 ]
COLD=wind[ wind$perc<.75 ]


cpgi = readBed("/data/akalin/Base/Annotation/hg19/CpGi/150114.UCSC.CpG.hg19.bed",track.line = 1)
cpgi.nonhot = cpgi[cpgi  %outside% wind.inter[["0.99-1"]] ]
regions = wind.inter
regions[["nonHOT-CpGi"]] = cpgi.nonhot
regions[["HOT-CpGi"]] = subsetByOverlaps(cpgi, wind.inter[["0.99-1"]])
regions[["HOT-nonCpGi"]] = wind.inter[["0.99-1"]][wind.inter[["0.99-1"]]  %outside% cpgi ]
regions[["CpGi"]] = cpgi


#############################################################################
## DNAseq-seq peaks
#############################################################################


dnaseq_files = "/data/akalin/Projects/AAkalin_HOTRegions/Results/Reviews_NAR/files_dnaseq_K562.txt"
dnaseq_meta = read.table("/data/akalin/Projects/AAkalin_HOTRegions/Results/Reviews_NAR/metadata_dnaseq_K562.tsv",
                         sep="\t", header = TRUE, stringsAsFactors = FALSE)
dnaseq_meta = dnaseq_meta[which(dnaseq_meta$File.Status=="released"),]
dnaseq_meta$File.accession

dir="/data/akalin/Projects/AAkalin_HOTRegions/Results/Reviews_NAR/"
bedfiles = paste0(dir, dnaseq_meta$File.accession, ".bed.gz")
library('genomation')
peaks = lapply(bedfiles, readNarrowPeak)
peaks = lapply(peaks, function(x) keepStandardChromosomes(x, "Homo_sapiens", pruning.mode="tidy"))
dnasepeaks = lapply(peaks, function(x) sort(x, by = ~signalValue, decreasing = TRUE))


dnasepeaks = lapply(dnasepeaks, function(x){
  x$rank = 1:length(x)
  return(x)
})

dnaseq.on.regions=lapply(dnasepeaks, function(x){
  return(list(
    HOT=subsetByOverlaps(x, HOT),
    MILD=subsetByOverlaps(x, MILD),
    COLD=subsetByOverlaps(x, COLD)
    ))
})

ranks.all.quantile = quantile(dnasepeaks[[1]]$rank, probs=seq(0, 1, .05))

peaks.ranks.all.quantile = lapply(2:length(ranks.all.quantile), function(i){
  prev = ranks.all.quantile[i-1]
  current = ranks.all.quantile[i]
  dnasepeaks[[1]][which(dnasepeaks[[1]]$rank <= current & dnasepeaks[[1]]$rank > prev), ]
})
names(peaks.ranks.all.quantile) = names(ranks.all.quantile)[2:length(ranks.all.quantile)]

peaks.ranks.all.quantile.len = sapply(peaks.ranks.all.quantile, length)


peaks.ranks.all.quantile.HOT = lapply(peaks.ranks.all.quantile,
                                      function(x) subsetByOverlaps(x, HOT)
)
peaks.ranks.all.quantile.HOT.len = sapply(peaks.ranks.all.quantile.HOT, length)

peaks.ranks.all.quantile.MILD = lapply(peaks.ranks.all.quantile,
                                      function(x) subsetByOverlaps(x, MILD)
)
peaks.ranks.all.quantile.MILD.len = sapply(peaks.ranks.all.quantile.MILD, length)

peaks.ranks.all.quantile.COLD = lapply(peaks.ranks.all.quantile,
                                      function(x) subsetByOverlaps(x, COLD)
)
peaks.ranks.all.quantile.COLD.len = sapply(peaks.ranks.all.quantile.COLD, length)


output.matrix1 = matrix(0, 3, length(peaks.ranks.all.quantile), 
                       dimnames = list(c("HOT", "MILD", "COLD"),
                                       names(peaks.ranks.all.quantile)))
output.matrix1[1,] = peaks.ranks.all.quantile.HOT.len
output.matrix1[2,] = peaks.ranks.all.quantile.MILD.len
output.matrix1[3,] =peaks.ranks.all.quantile.COLD.len

tmmp=c('0-5',
'5-10',
'10-15',
'15-20',
'20-25',
'25-30',
'30-35',
'35-40',
'40-45',
'45-50',
'50-55',
'55-60',
'60-65',
'65-70',
'70-75',
'75-80',
'80-85',
'85-90',
'90-95',
'95-100')
colnames(output.matrix1) =tmmp 

pdf("~/projects/HOT-regions/Review_NAR/DNAseqHOTMILDCOLD.pdf", width=6, height = 5)
barplot(output.matrix1,
        ylab= "Number of regions on accessible sites",
        xlab = "\n\nTop most accessible sites [percentile]",
        col=c("#ff0000ff", "#00a900ff", "#000080ff"), 
        las=2)
legend("topright", 
       legend = c("HOT", "MILD", "COLD"), 
       fill = c("#ff0000ff", "#00a900ff", "#000080ff"))
dev.off()



