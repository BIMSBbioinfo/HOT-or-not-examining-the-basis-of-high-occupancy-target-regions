
options(scipen=999) #disable scientific notation
.libPaths(c("/home/kwreczy/Rlibs/3.4/"))

################
############ IMPORTS
###############

library(GenomicRanges)
library(genomation)
library(tools)
library(rtracklayer)
library(parallel)

################
############ FUNCTIONS
###############

setwd("~/projects/HOT-regions/Define_HOT_regions/")
source("~/projects/HOT-regions/Define_HOT_regions/peaker_adapted.R")

gr2df = function(gr){
  
  data.frame(seqnames=seqnames(gr),
             start=start(gr)-1,
             end=end(gr),
             names=c(rep(".", length(gr))),
             score=gr$scores,
             strand=strand(gr))
}

gr2df.narrowPeak = function(gr){
  data.frame(seqnames=as.character(seqnames(gr)),
             starts=start(gr)-1,
             ends=end(gr),
             names=c(rep(".", length(gr))),
             scores=c(rep("0", length(gr))),
             strands=c(rep("*", length(gr))),
             signalValue=c(rep("0", length(gr))),
             pValue=c(rep("0", length(gr))),
             qValue=c(rep("0", length(gr))),
             peak=gr$peak
  )
}

# copied pasted from genomation R package
readBigWig = function(target, windows=NULL, ...){
  if(class(windows) != 'GRanges')
    stop('windows argument needs to be a GRanges object')
  if(is.null(windows)){
    bw = import(target)
  }else{
    bw = import(target, which=windows)
  }
  if(length(bw) == 0)
    stop('There are no ranges selected')
  covs = coverage(bw, weight=bw$score)
  return(covs)
}

parallel_wigToBigWig = function(i, wigfiles=wigfiles, seqinfo=seqinfo){
  print(i)
  wigToBigWig(wigfiles[i], seqinfo)
}

bw.and.gff2Narrowpeak_df = function(bigwigfile, gfffile, add.chr=TRUE, seqinfo=NULL)
  ## add.chr - adds prefix "chr" to names of chromosomes
{
  windows = gffToGRanges(gfffile) 
  windows = sort(sortSeqlevels(windows))
  if(add.chr==TRUE){
    seqlevels(windows) = paste("chr",seqlevels(windows), sep="")
  }
  if(!is.null(seqinfo)){
    seqinfo(windows) <- seqinfoo[seqnames(seqinfo(windows))]
    windows = trim(windows)
  }
  
  cov = readBigWig(bigwigfile,windows)
  
  summits = c()
  chrms = seqlevels(windows)
  # for each chromosome
  for(chr in chrms){
    wind = windows[seqnames(windows)==chr,]
    maxs=viewWhichMaxs(Views(cov[[chr]],start(wind), end(wind)), na.rm=TRUE) #summits within windows
    summits = c(summits, maxs)
  }
  windows$peak = end(windows)-summits
  df = gr2df.narrowPeak(windows)  
  df
}

################
############ MAIN
###############

###### mouse (mm9) 

meta.mm9 = read.table("/data/akalin/Projects/AAkalin_HOTRegions/Data/mouse_TF_ENCODE/NarrowPeak/metadata.tsv",
                      header=TRUE,sep="\t")
files = paste0(meta.mm9$File.accession, ".bed.gz")

# convert all TFs fom ENCODE in NarrowPeak format to GRangesList
tf.path = "/data/akalin/Projects/AAkalin_HOTRegions/Data/mouse_TF_ENCODE/NarrowPeak/"
tf_list_of_gr <- sapply(paste0(tf.path, files), function(x) readNarrowPeak(x))
tf_grl <- GRangesList(tf_list_of_gr)
names(tf_grl) <- meta.mm9$File.accession

pcov500 = getPeakSummitcov(tf_grl)
peakAnchorsW05kFl=getPeakRle(pcov500,filterNearBy=1000)
write.table(gr2df(peakAnchorsW05kFl), file="/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/HOTmm9.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)

#' Remove polymerase binding sites to call HOT regions.
mm9.polymerase = c("POLR2AphosphoS2-mouse", "POLR2A-mouse")
tf_grl.nopol2a = tf_grl[-which( meta.mm9$Experiment.target %in% mm9.polymerase)]
pcov500 = getPeakSummitcov(tf_grl.nopol2a)
peakAnchorsW05kFl=getPeakRle(pcov500,filterNearBy=1000)
peakAnchorsW05kFl =  peakAnchorsW05kFl[order(-peakAnchorsW05kFl$scores)]




pdf("/home/kwreczy/projects/HOT-regions/Define_HOT_regions/hist_chipseqcount_mm9.pdf", width = 7, height = 5)
hot = peakAnchorsW05kFl
mcols(hot )$perc=ecdf(hot $scores)(hot $scores)
hhot = hot [ mcols(hot )$perc>=.99 ]

hist(countOverlaps(hot,tf_grl.nopol2a ,ignore.strand=TRUE),
     breaks=30,
     col="blue",
     main="Histogram of ChIP-seq peak counts per region",
     xlab="Number of ChIP-seq TF peaks")
abline(v=quantile(mcols(hot )$scores,probs = c(.75)),
       col="red")
text(8,100000, "75th percentile", col = "red",  srt=90, pos=3)
abline(v=quantile(mcols(hot )$scores,probs = c(.99)),
       col="red")
text(28,100000, "99th percentile", col = "red",  srt=90, pos=3)
dev.off()


###### human (hg19)

# HOT regions based on only peaks from Uniform UCSC track.
library(readr)
library(GenomicRanges)
bednarropeak2gr <- function(file){
  df <- read_delim(file, col_types="ciicicnnni", delim="\t",col_names = FALSE,
                   locale=locale(decimal_mark = ",", grouping_mark = "_"))
  colnames(df) <- c("chrom","start", "end", "name", "score", "strand","signalValue","qValue","pValue","peak")
  g = makeGRangesFromDataFrame(
    df,
    keep.extra.columns=FALSE,
    starts.in.df.are.0based=FALSE,
    ignore.strand=TRUE)
  g$name = df$name
  g$score =df$score
  g$signalValue = df$signalValue
  g$qValue = df$qValue
  g$pValue = df$pValue
  g$peak = df$peak
  g
}
# convert all TFs in NarrowPeak format to GRangesList
meta = read.table("/data/akalin/Projects/AAkalin_HOTRegions/Data/wgEncodeAwgTfbsUniform/UCSC_Uniform_TFBS_metafile.txt",
                  header=T, sep="\t")
tf.path = "/data/akalin/Projects/AAkalin_HOTRegions/Data/wgEncodeAwgTfbsUniform/narrowPeak/"
l = paste0(tf.path, paste0(meta$tableName, ".narrowPeak.gz"))
tf_list_of_gr <- mclapply(l, bednarropeak2gr, mc.cores=20)
tf_grl <- GRangesList(tf_list_of_gr)
names(tf_grl) <- meta$Factor
#saveRDS(tf_grl, "/data/akalin/kasia/projects/HOT_DATA/rds/hg19_tfbs_uniform.rds")
#tf_grl <- readRDS("/data/akalin/kasia/projects/HOT_DATA/rds/hg19_tfbs_uniform.rds")

# remove pol2a
polymerase = c("POLR3G","POLR2A")
meta.pol2a.inx = which(meta$Factor %in% polymerase)
tf_grl.nopol2a = tf_grl[-meta.pol2a.inx]

pcov500 = getPeakSummitcov(tf_grl.nopol2a, p=501)
peakAnchorsW05kFl=getPeakRle(pcov500,filterNearBy=2000)

write.table(gr2df(peakAnchorsW05kFl[order(-peakAnchorsW05kFl$scores)]), 
            file="/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/HOThg19_filterNearBy2kb_Uniform_nopol2a_summits.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)


tmp  = findOverlaps(hot, tf_grl.nopol2a)
tmp1 = as.data.table(tmp)      
tmp2 = tmp1[,.(COUNT=.N), by=.(queryHits)]
tmp3 = cbind(tmp2, score=hot$scores[tmp2$queryHits])
       


pdf("/home/kwreczy/projects/HOT-regions/Define_HOT_regions/hist_chipseqcount_hg19.pdf", width = 7, height = 5)
hot = peakAnchorsW05kFl
mcols(hot )$perc=ecdf(hot $scores)(hot $scores)
hhot = hot [ mcols(hot )$perc>=.99 ]

hist(countOverlaps(hot,tf_grl.nopol2a ,ignore.strand=TRUE),
     breaks=30,
     col="blue",
     main="Histogram of ChIP-seq peak counts per region",
     xlab="Number of ChIP-seq TF peaks")
abline(v=quantile(mcols(hot )$scores,probs = c(.75)),
       col="red")
text(33,150000, "75th percentile", col = "red",  srt=90, pos=3)
abline(v=quantile(mcols(hot )$scores,probs = c(.99)),
       col="red")
text(200,150000, "99th percentile", col = "red",  srt=90, pos=3)
dev.off()



###### c.elegans (ce10)

## Seqinfo, info about length of chromosomes
# (note: before ce6.)
ce10.chrom.sizes = read.table("http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.chrom.sizes",
                              stringsAsFactors=FALSE, colClasses=c("character","integer"),
                              col.names=c("chrom","size"))
g = c("RHet", "LHet", "R", "L")
chrms = as.vector(sapply(g, function(y) sapply(1:5, function(x) paste0("chr", x, y))))
ce10.chrom.sizes = rbind(ce10.chrom.sizes,data.frame(chrom=c("chrU","chrUextra", chrms), 
                                                     size=1e8))
# 1e8 is just a random big number. it shouldn't have any impact, because
# anyways regions will have to be trimmed to chromosomes sizes.
# without this seqinfo object I am not able to use wigToBigWig function.
ce10.chrom.sizes[which(ce10.chrom.sizes$chrom=="chrV"),]$size <- 1e8 
ce10.chrom.sizes[which(ce10.chrom.sizes$chrom=="chrII"),]$size <- 1e8
ce10.chrom.sizes[which(ce10.chrom.sizes$chrom=="chrIII"),]$size <- 1e8
ce10.chrom.sizes[which(ce10.chrom.sizes$chrom=="chrX"),]$size <- 1e8
ce10.chrom.sizes = rbind(ce10.chrom.sizes, data.frame(chrom="chr1", size=ce10.chrom.sizes[ce10.chrom.sizes$chrom=="chrI",]$size))
ce10.chrom.sizes = rbind(ce10.chrom.sizes, data.frame(chrom="chr2", size=ce10.chrom.sizes[ce10.chrom.sizes$chrom=="chrII",]$size))
ce10.chrom.sizes = rbind(ce10.chrom.sizes, data.frame(chrom="chr3", size=ce10.chrom.sizes[ce10.chrom.sizes$chrom=="chrIII",]$size))
ce10.chrom.sizes = rbind(ce10.chrom.sizes, data.frame(chrom="chr4", size=ce10.chrom.sizes[ce10.chrom.sizes$chrom=="chrIV",]$size))
ce10.chrom.sizes = rbind(ce10.chrom.sizes, data.frame(chrom="chr5", size=ce10.chrom.sizes[ce10.chrom.sizes$chrom=="chrV",]$size))

#e.g. 3 instead of chr3
without.chr.chrom = substring(ce10.chrom.sizes$chrom, 4, nchar(ce10.chrom.sizes$chrom)) 
temp = data.frame(without.chr.chrom, ce10.chrom.sizes$size)
colnames(temp) <- c("chrom","size")
ce10.chrom.sizes <- rbind(ce10.chrom.sizes, temp)
# seqinfo object
seqinfoo = seqinfo(ce10.chrom.sizes$chrom,
                   seqlengths=ce10.chrom.sizes$size+100,
                   isCircular=NA, genome="ce10")

## Calculate point source values (pick summits)
## by using .wig and .gff3 files and based on them generate .narrowpeak
# convert wig to bigWig
#meta = read.table("/data/akalin/Projects/AAkalin_HOTRegions/Data/c_elegans_TF_ENCODE/metadata.tsv", sep="\t", header=T)
tf.path = "/data/akalin/Projects/AAkalin_HOTRegions/Data/c_elegans_TF_ENCODE/"
l <- list.files(path = tf.path, pattern = "*.wig")
wigfiles <- sapply(l, function(x) paste(tf.path, x,sep=""))

#mclapply(1:length(wigfiles), parallel_wigToBigWig, wigfiles=wigfiles,seqinfo=seqinfoo, mc.cores = 20)#why it doesnt work?
for(i in 1:length(wigfiles)){
  a = try(wigToBigWig(wigfiles[i], seqinfoo))
  if(is(a,"try-error")) {
    print(a)
    next 
  }
}
tf.path = "/data/akalin/Projects/AAkalin_HOTRegions/Data/c_elegans_TF_ENCODE_TF_ENCODE/"
l <- list.files(path = tf.path, pattern = "*.wig")
wigfiles <- sapply(l, function(x) paste(tf.path, x,sep=""))
#mclapply(1:length(wigfiles), parallel_wigToBigWig, wigfiles=wigfiles,seqinfo=seqinfoo, mc.cores = 20)
tf.path = "/data/akalin/Projects/AAkalin_HOTRegions/Data/c_elegans_TF_ENCODE/"
l <- list.files(path = tf.path, pattern = "*combined.bw")
bwfiles <- sapply(l, function(x) paste(tf.path, x,sep=""))

# take only ip and dont take input. take files that dont have "input" in their file names
bwfiles1 = c()
for(i in 1:length(bwfiles)){
  if(!grepl("nput",bwfiles[i])){
    bwfiles1 = c(bwfiles1,bwfiles[i] )
  }
}
bwfiles = bwfiles1
bwfilenames <- basename(file_path_sans_ext(bwfiles))
l <- list.files(path = tf.path, pattern = "*.gz")
gfffiles <- sapply(l, function(x) paste(tf.path, x,sep=""))
gfffilenames <- basename(file_path_sans_ext(file_path_sans_ext(gfffiles)))

# for some .bw there is no .gff
bwgff = data.frame(stringsAsFactors = FALSE)
for(i in 1:length(bwfilenames)){
  a = which(grepl(bwfilenames[i], gfffilenames))
  if(length(a)>=1){
    bwgff <- rbind(bwgff, data.frame(bwfiles[i], gfffiles[a], stringsAsFactors = FALSE))
  }
}
bwfilesnames <- basename(file_path_sans_ext(file_path_sans_ext(file_path_sans_ext(bwgff[,1]))))
narrowpeakfiles = paste(tf.path,bwfilesnames,".narrowPeak", sep="")
for(i in 1:nrow(bwgff)){
  print(i)
  df.narrow = try(bw.and.gff2Narrowpeak_df(as.character(bwgff[i,][1]), as.character(bwgff[i,][2]), add.chr=FALSE))
  print(unique(df.narrow[,1]))
  if(is(df.narrow,"try-error")) {
    print(i)
    df.narrow = try(bw.and.gff2Narrowpeak_df(as.character(bwgff[i,][1]), as.character(bwgff[i,][2]), add.chr=TRUE))
    if(is(df.narrow,"try-error")) {
      print(df.narrow)
    }
  }  
}
l <- list.files(path = tf.path, pattern = "*.narrowPeak")
l <- sapply(l, function(x) paste(tf.path, x,sep=""))
tf_list_of_gr <- sapply(l, function(x) readNarrowPeak(x))

tf_grl <- GRangesList(tf_list_of_gr)
pcov500 = getPeakSummitcov(tf_grl)
peakAnchorsW05kFl=getPeakRle(pcov500,filterNearBy=1000)
#saveRDS(peakAnchorsW05kFl,"~/HOT/HOT_DATA/rds/peakAnchorsW05kFl.rds")
write.table(gr2df(peakAnchorsW05kFl), file="/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTce10.bed", quote=F, sep="\t", row.names=F, col.names=F)


###### d.melanogaster (dm3)
#meta = read.table("/data/akalin/Projects/AAkalin_HOTRegions/Data/d_melanogaster_TF_ENCODE/metadata.tsv", sep="\t", header=T)
# Seqinfo, info about length of chromosomes
dm3.chrom.sizes = read.table("http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.chrom.sizes",
                             stringsAsFactors=FALSE, colClasses=c("character","integer"),
                             col.names=c("chrom","size"))

# Note: Ensembl uses chromosome names of "1" where UCSC uses: "chr1"
without.chr.chrom = substring(dm3.chrom.sizes$chrom, 4, nchar(dm3.chrom.sizes$chrom)) #e.g. 3 instead of chr3

# some chromosomes are divided into two parts, like chr3 into chr3R and chr3L
chr3.size = dm3.chrom.sizes[dm3.chrom.sizes$chrom=="chr3L",]$size + dm3.chrom.sizes[dm3.chrom.sizes$chrom=="chr3R",]$size

# seqinfo object
seqinfoo = Seqinfo(c(dm3.chrom.sizes$chrom,
                     without.chr.chrom,
                     "chrdmel_mitochondrion_genome", 
                     "dmel_mitochondrion_genome", 
                     "3", 
                     "chr3"),
                   seqlengths=c(dm3.chrom.sizes$size+140,
                                dm3.chrom.sizes$size+140,
                                19656,
                                19656,
                                chr3.size,
                                chr3.size),
                   isCircular=NA, genome="dm3")
real_seqinfo = Seqinfo(c(dm3.chrom.sizes$chrom,
                         without.chr.chrom,
                         "chrdmel_mitochondrion_genome", 
                         "dmel_mitochondrion_genome", 
                         "3", 
                         "chr3"),
                       seqlengths=c(dm3.chrom.sizes$size,
                                    dm3.chrom.sizes$size,
                                    19517,
                                    19517,
                                    chr3.size,
                                    chr3.size),
                       isCircular=NA, genome="dm3")


## Calculate point source values like from NarrowPeak format (pick summits)
## by using .wig and .gff3 files
# peaks in .gff3 files are calculated by MACS
# convert wig to bigWig
tf.path = "/data/akalin/Projects/AAkalin_HOTRegions/Data/d_melanogaster_TF_ENCODE/"
l <- list.files(path = tf.path, pattern = "*.wig")
wigfiles <- sapply(l, function(x) paste(tf.path, x,sep=""))
#mclapply(1:length(wigfiles), parallel_wigToBigWig, wigfiles=wigfiles,seqinfo=seqinfoo, mc.cores = 20)
tf.path = "/data/akalin/Projects/AAkalin_HOTRegions/Data/d_melanogaster_TF_ENCODE/"
l <- list.files(path = tf.path, pattern = "*.bw")
bwfiles <- sapply(l, function(x) paste(tf.path, x,sep=""))
bwfilenames <- basename(file_path_sans_ext(bwfiles))
bwfilenames.short <- sapply(bwfilenames, function(x) substring(x, 1, nchar(x)-8)) #wihtout _density at the end of the file
bwfilenames.shorter <- substring(bwfilenames.short, 1,7)

tf.path = "/data/akalin/Projects/AAkalin_HOTRegions/Data/d_melanogaster_TF_ENCODE/"
allfiles <- list.files(path = tf.path, pattern = "*")

gfffiles.basename = unlist(sapply(bwfilenames.shorter, function(x) allfiles[which(grepl(paste(".*",x,".*.gz$", sep=""),allfiles,fixed = FALSE))]))
gfffiles <- paste(tf.path,gfffiles.basename,sep="")
# check if all gff files exist
length(gfffiles) == length(which(!(sapply(gfffiles, function(x) file.exists(paste(tf.path,x,sep=""))))))
narrowpeakfiles = paste(tf.path,bwfilenames.short,".narrowPeak", sep="")
for(i in 1:length(bwfiles)){
  print(i)
  df.narrow = try(bw.and.gff2Narrowpeak_df(bwfiles[i], gfffiles[i], add.chr=TRUE, seqinfo=real_seqinfo))
  if(class(df.narrow)!="try-error"){
    df.narrow$seqnames = as.character(df.narrow$seqnames)
    m = df.narrow[df.narrow$seqnames=="chrdmel_mitochondrion_genome",]$seqnames
    if(length(m)>0){
      df.narrow[df.narrow$seqnames=="chrdmel_mitochondrion_genome",]$seqnames = "chrM"
    }
    print(unique(df.narrow[,1]))
  }
  
  if(class(df.narrow)=="try-error") {  
    print(i)
    df.narrow = try(bw.and.gff2Narrowpeak_df(bwfiles[i], gfffiles[i], add.chr=FALSE,seqinfo=real_seqinfo))
    print("add chr false")
    print(unique(df.narrow[,1]))
    df.narrow$seqnames = as.character(df.narrow$seqnames)
    df.narrow$seqnames = paste0("chr",df.narrow$seqnames)
    m = df.narrow[df.narrow$seqnames=="chrdmel_mitochondrion_genome",]$seqnames
    if(length(m)>0){
      df.narrow[df.narrow$seqnames=="chrdmel_mitochondrion_genome",]$seqnames = "chrM"
    }
    print(unique(df.narrow[,1]))
    if(is(df.narrow,"try-error")) {
      print(df.narrow)
    }
  }  
  write.table(df.narrow, narrowpeakfiles[i], col.names = F, row.names = F, quote = F, sep="\t")
}

# convert all TFs in NarrowPeak format to GRangesList
l <- list.files(path = tf.path, pattern = "*.narrowPeak")
l <- sapply(l, function(x) paste(tf.path, x,sep=""))
tf_list_of_gr <- sapply(l, function(x) readNarrowPeak(x))

tf_grl <- GRangesList(tf_list_of_gr)
pcov500 = getPeakSummitcov(tf_grl)
peakAnchorsW05kFl=getPeakRle(pcov500,filterNearBy=1000)
write.table(gr2df(peakAnchorsW05kFl), file="/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTdm3.bed", quote=F, sep="\t", row.names=F, col.names=F)


