#'---
#'title: "Rloops/G4/Methylation vs HOT score"
#'author: "Katarzyna Wreczycka"
#'date: "`r format(Sys.time(), '%B %d, %Y')`"
#'output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    fig_width: 7
#'    fig_height: 7
#'---
#'
#'

library(genomation)
library(Rsamtools)
library(ggplot2)
library("gridExtra")
library(reshape)
library(cowplot)
library(plyr)
library(reshape)


my_theme = theme(
  strip.text.x = element_text(size=8),
  axis.text.x = element_text(angle = 45, hjust = 1, size=8),
  axis.text.y = element_text(hjust = 1, size=8),
  axis.title=element_text(size=8,face="plain"),
  axis.line = element_line(size = 1, linetype = "solid", colour = "black"), 
  panel.grid.major = element_line(colour = "#E6E6E6"))


################
############ FUNCTIONS
###############

source("~/projects/HOT-regions/basic_functions.R")
#' Blacklisted regions for hg19
black = readBed("/home/buyar/datasets/genomes/hg19/wgEncodeDacMapabilityConsensusExcludable.bed")
BS.genomeSpecies=get.BS.genomeSpecies("hg19")
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

wins.vec = function(vec, rng=c(0, 100)){
  hi.th=quantile(vec,rng[2]/100,na.rm=TRUE)
  vec[vec>=hi.th]=hi.th
  lo.th=quantile(vec,rng[1]/100,na.rm=TRUE)
  vec[vec<=lo.th]=lo.th
  vec
}
#Emulate ggplot2 default color palette
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color1 = ggplotColours(n=3)[3] # blue color like blue color in ggplot pallette


# function for IP/control enrichment over pre-defined regions
enrich<-function(signal,background,regions){
  require(Rsamtools)
  param <- ScanBamParam(which=regions)
  cnts.sg=countBam(signal, param=param)$records
  cnts.bg=countBam(background, param=param)$records
  
  message("normalizing... ",date())
  command=paste0("/home/aakalin/.guix-profile/bin/samtools idxstats ",signal," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
  sg.mapped=as.numeric(try(system(command,intern=TRUE)))
  command2=paste0("/home/aakalin/.guix-profile/bin/samtools idxstats ",background," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
  bg.mapped=as.numeric(try(system(command2,intern=TRUE)))
  
  cnts.sg=1e6*cnts.sg/sg.mapped
  cnts.bg=1e6*cnts.bg/bg.mapped
  log2((cnts.sg+1)/(cnts.bg+1))
}

calc.signal<-function(signal,regions){
  require(Rsamtools)
  param <- ScanBamParam(which=regions)
  cnts.sg=Rsamtools::countBam(signal, param=param)$records
  
  message("normalizing... ",date())
  command=paste0("/home/aakalin/.guix-profile/bin/samtools idxstats ",signal," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
  sg.mapped=as.numeric(try(system(command,intern=TRUE)))
  
  cnts.sg=1e6*cnts.sg/sg.mapped
  log2(cnts.sg+1)
}

.jets<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

################
############ READ DATA
###############

#' Read HOT summits in hg19
#' 
source("~/projects/HOT-regions/basic_functions.R")
BS.genomeSpecies=get.BS.genomeSpecies("hg19")
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

hg19 = processHOTsummits("/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/HOThg19_filterNearBy2kb_Uniform_nopol2a_summits.bed",
                         2000, BS.genomeSpecies)
# remove unnecessary chromosomes
hg19 = keepSeqlevels(hg19,
                     unique(as.character(seqnames(hg19))) )
black = readBed("/home/buyar/datasets/genomes/hg19/wgEncodeDacMapabilityConsensusExcludable.bed")
hg19noblack = hg19[which(hg19 %outside% black)]
hg19noblack = hg19noblack[order(-hg19noblack$perc)]
hot=hg19noblack
wind = hg19noblack
wind = add.intervals(wind)
wind = keepSeqlevels(wind, unique(as.character(seqnames(wind))) )


################
############ MAIN
###############

#" ## Nadel et al 2015 GSE68948 and Lim et al 2015 GSE57353
#' Read paths to drip-seq samples
#' p=read.table("~/projects/HOT-regions/DRIP-seq/GSE68948.txt", sep="\t", header=TRUE, stringsAsFactors =FALSE)
#' s=read.table("~/projects/HOT-regions/DRIP-seq/GSE57353.txt", sep="\t", header=TRUE, stringsAsFactors =FALSE)
#' inp = rbind(p, s)
#' print(inp)
#' 
#' #' Drip-seq peaks attached to the publications
#' GSE57353_controls_path = "/home/kwreczy/projects/HOT-regions/DRIP-seq/GSE57353_controls_drip.bed"
#' GSE68948_HEKminus10_path = "/home/kwreczy/projects/HOT-regions/DRIP-seq/GSE68948_HEKminus10.txt"
#' GSE68948_IMR90minus10_path ="/home/kwreczy/projects/HOT-regions/DRIP-seq/GSE68948_IMR90minus10.txt"
#' 
#' GSE70189 = paste0("/data/akalin/kasia/GSE70189",c("GSM1720615_NT2_DRIP_1_peaks.bed.gz",
#'                                                   "GSM1720616_NT2_DRIP_2_peaks.bed.gz",
#'                                                   "GSM1720617_NT2_DRIP_RNaseA_peaks.bed.gz",
#'                                                   "GSM1720618_NT2_DRIP_RNaseH_peaks.bed.gz"	,
#'                                                   "GSM1720619_K562_DRIP_peaks.bed.gz"	
#' ))
#' signal = c("IMR90_RDIP", "HEK293T_RDIP","GSE57353Control1", "GSE57353Control3")
#' background = c("IMR90_Input","HEK293T_Input","GSE57353Digestedinput","GSE57353Digestedinput")
#' ip.bams = sapply(signal, 
#'                  function(x) 
#'                    paste0(hdir, inp[which(inp$name==x),]$dirname, "/MAPPED/Bowtie/", inp[which(inp$name==x),]$run, ".sorted.bam"))
#' ctrl.bams = sapply(background, 
#'                    function(x) 
#'                      paste0(hdir, inp[which(inp$name==x),]$dirname, "/MAPPED/Bowtie/", inp[which(inp$name==x),]$run, ".sorted.bam"))
#result=mcmapply(enrich,ip.bams,ctrl.bams,
#                   MoreArgs = list(regions = hot),
#                   mc.cores=4) 
#colnames(result)=names(ip.bams)
#result = cbind(as.data.frame(result), perc=hot$perc)
# hot = add.intervals(hot)
# #result = cbind(as.data.frame(result), x=hot$x)
#save.logFC.rds="/data/akalin/kwreczy/projects/HOT_DATA/rds/log2Enrichment.dripseq_GSE68948_GSE57353.on.HOT_Uniform_nopol2a.rds"
# #saveRDS(result,
# #        save.logFC.rds)

save.logFC.rds="/data/akalin/kwreczy/projects/HOT_DATA/rds/log2Enrichment.dripseq_GSE68948_GSE57353.on.HOT_Uniform_nopol2a.rds"
result = readRDS(save.logFC.rds)
res = result[,c(1,2,3,4,6)]
colnames(res) <- c("IMR90 cell line","HEK293T cell line","Primary\nfibroblasts rep1","Primary\nfibroblasts rep2", "x")
result.GSE68948_GSE57353 = res
dfmelt_GSE68948_GSE57353<-melt(res, measure.vars = 1:4)  



# GSE70189.table = read.table("/data/akalin/kwreczy/GSE70189/GSE70189.txt", sep="\t", header=TRUE, stringsAsFactors =FALSE)
# GSE70189.table.hg19 = GSE70189.table[which(GSE70189.table$scientific_name=="Homo sapiens"),]
# ip  = c("SRR2075681","SRR2075682") #biological replicates
# ctrl  = c("SRR2075684")
# GSE70189.ip.bams = paste0(paste0("/data/akalin/kwreczy/GSE70189/",ip), ".sorted.bam")
# GSE70189.ctr.bams = paste0(paste0("/data/akalin/kwreczy/GSE70189/",c(ctrl,ctrl)), ".sorted.bam")
# result=mcmapply(enrich,
#                 GSE70189.ip.bams,
#                 GSE70189.ctr.bams,
#                    MoreArgs = list(regions = wind),
#                    mc.cores=2)
# colnames(result) <- GSE70189.ip.bams
# result <- cbind(result, perc=wind$perc, x=wind$x)
# result = as.data.frame(result,stringsAsFactors = FALSE)
# result$x=NULL
# result$x=as.factor(as.character(wind$x))
#  saveRDS(result,
#         "/data/akalin/kwreczy/projects/HOT_DATA/rds/log2signal.dripseqGSE70189.ctrlrnaseh.rds")
result=readRDS("/data/akalin/kwreczy/projects/HOT_DATA/rds/log2signal.dripseqGSE70189.ctrlrnaseh.rds")
result.GSE70189 = result
colnames(result.GSE70189) <- c("NT2 cell line rep1" ,   "NT2 cell line rep2","perc", "x")
result.GSE70189 = result.GSE70189[,c(1,2,4)]
dfmelt_GSE70189<-melt(result.GSE70189, measure.vars = 1:2)  

rangee = c(min(dfmelt_GSE68948_GSE57353$value,dfmelt_GSE70189$value), max(dfmelt_GSE68948_GSE57353$value,dfmelt_GSE70189$value))

dripseq_GSE68948_GSE57353<-ggplot(dfmelt_GSE68948_GSE57353, aes(x=x, y=value, fill=variable))+
  #geom_boxplot(colour="black", fill="steelblue")+
  #stat_boxplot(geom ='errorbar') + 
  geom_boxplot(colour="black", fill="steelblue", 
               outlier.size = 0.5)+
  ylim(rangee)+
  facet_grid(.~variable)+
  labs(y="DRIP-seq log2(IP/Input)\n", x="\nTF occupancy percentiles") +
  #coord_cartesian(ylim = rangee)+
  my_theme 

dripseq_GSE70189.1<-ggplot(dfmelt_GSE70189, aes(x=x, y=value, fill=variable))+
  #stat_boxplot(geom ='errorbar') + 
  geom_boxplot(colour="black", fill="steelblue", 
               outlier.size = 0.5)+
  ylim(rangee)+
  #geom_boxplot(colour="black", fill="steelblue")+
  facet_grid(.~variable)+
  labs(y="DRIP-seq log2(IP/(IP+RNaseH))\n", x="\nTF occupancy percentiles") +
  my_theme 


#### ce10 DRIP-seq

result = readRDS("/data/akalin/kwreczy/projects/HOT_DATA/rds/log2Enrichment.dripseq_ce10_swapped.rds")
dfmelt_input<-melt(result[,c(1,4)], measure.vars = 1)  
dripseq_input_ce10<-ggplot(dfmelt_input, aes(x=x, y=value))+
  #stat_boxplot(geom ='errorbar') + 
  geom_boxplot(colour="black", fill="steelblue", 
               outlier.size = 0.5)+
  #ylim(-2.1, 2.1)+
  #geom_boxplot(colour="black", fill="steelblue")+
  labs(y="DRIP-seq log2(IP/Input)\n", x="\nTF occupancy percentiles") +
  my_theme 

#### yeast DRIP-seq

calc.signal<-function(signal,regions){
  require(Rsamtools)
  param <- ScanBamParam(which=regions)
  cnts.sg=Rsamtools::countBam(signal, param=param)$records
  
  message("normalizing... ",date())
  command=paste0("/home/aakalin/.guix-profile/bin/samtools idxstats ",signal," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
  sg.mapped=as.numeric(try(system(command,intern=TRUE)))
  
  cnts.sg=1e6*cnts.sg/sg.mapped
  log2(cnts.sg+1)
}

Integrate_Yeast_Hyper_Drip = function(){
  
  library(GenomicRanges)
  library(genomation)
  library(rtracklayer)
  hyper_path = list.files("/data/akalin/Projects/AAkalin_HOTRegions/AccessoryData/Teytelman_2013_PNAS_HyperRegs", full.names=TRUE, pattern='bed')
  hyper = genomation::readGeneric(hyper_path)
  
  drip_path = list.files("/data/akalin/Projects/AAkalin_HOTRegions/AccessoryData/SRR3504389_Wahba_2016_GenDev_SCev.Dripseq"
                         , full.names=TRUE,
                         pattern='bam$', recursive=TRUE)
  
  yeast_annot = file.path("/data/akalin/Base",'Annotation','sacCer3','170208_sacCer3_UCSC_EnsemblGenes.gtf')
  #annot = RCAS::importGtf(yeast_annot, keepStandardChr=FALSE)
  annot = rtracklayer::import.gff(yeast_annot)
  annot = subset(annot, type=='exon')
  annot = annot[width(annot) > 10]
  sml = ScoreMatrixBin(drip_path, annot, bin.num=1, type='bam', strand.aware=FALSE)
  dat = data.frame(drip=as.vector(sml),
                   hyper = ifelse(countOverlaps(annot, hyper)>0,'Hyper','Not Hyper'))
  dat$hyper = factor(dat$hyper, levels=c('Not Hyper', 'Hyper'))

  #calc.signal(drip_path, annot) #yeast_dat_calc_signal #yeast_dat_calc_signal
  
}
#yeast_dat <- Integrate_Yeast_Hyper_Drip()
#saveRDS(yeast_dat, "/data/akalin/kwreczy/projects/HOT_DATA/rds/yeast_dat.rds")
yeast_dat <- readRDS("/data/akalin/kwreczy/projects/HOT_DATA/rds/yeast_dat.rds")
library(plyr)
yeast_dat$hyper = revalue(yeast_dat$hyper , c("Not Hyper"="Non-hyper", "Hyper"="Hyper"))


dfmelt.yeast_dat<-melt(yeast_dat, measure.vars = 1)  
dfmelt.yeast_dat$x = dfmelt.yeast_dat$hyper
yeast_dat_m<-ggplot(dfmelt.yeast_dat, aes(x=x, y=value))+
  #stat_boxplot(geom ='errorbar') + 
  geom_boxplot(colour="black", fill="steelblue", 
               outlier.size = 0.5)+
  ylim(0, 150)+
  labs(y="Average number of\nreads from DRIP-seq\nper base", x="\nOverlap with\nhyper-ChIPable\nregions") +
  my_theme 
  

#### hg19 G-quadr


log2enrich.G4 = readRDS("/home/kwreczy/projects/HOT-regions/Gquadruplex/log2enrich.G4.rds")


library(ggplot2)
library(reshape)
library(ggplot2)
colnames(log2enrich.G4) <- c("HaCaT cell line","NHEK cell line\nrep1", "NHEK cell line\nrep2","x")
dfmeltm<-melt(log2enrich.G4, measure.vars = 1:3)  
g4chipseq_m<-ggplot(dfmeltm, aes(x=x, y=value, fill=variable))+
  #stat_boxplot(geom ='errorbar') + 
  geom_boxplot(colour="black", fill="steelblue", 
               outlier.size = 0.5)+
  #ylim(-1, 1)+
  #geom_boxplot(colour="black", fill="steelblue")+
  facet_grid(.~variable)+
  labs(y="G4-ChIP-seq log2(IP/Input)\n", x="\nTF occupancy percentiles") +
  my_theme 

#### hg19 methylation

DNAme_WGBS.path = "/data/akalin/Base/RoadmapEpigenomics/Experiment/DNAme_WGBS/"
tabl = read.table(paste0(DNAme_WGBS.path,"EG.mnemonics.name.txt"), stringsAsFactors =FALSE)
# print(tabl)

FractionalMethylation.path = "/data/akalin/Base/RoadmapEpigenomics/Experiment/DNAme_WGBS/FractionalMethylation_bigwig/"
meth.files = paste0(FractionalMethylation.path, tabl[,1], "_WGBS_FractionalMethylation.bigwig" )


#' HOT regions
source("~/projects/HOT-regions/basic_functions.R")
#BS.genomeSpecies=get.BS.genomeSpecies("hg19")
require(rtracklayer)
hg19.seqinfo = SeqinfoForUCSCGenome("hg19")
processHOTsummits=function(fname,extend=1000,bs,add.chr=FALSE){
  require(genomation)
  require(GenomicRanges)
  require(rtracklayer)
  # read file
  hot=readBed(fname)
  if(add.chr)
    seqlevels(hot)=paste0("chr",seqlevels(hot))
  hot=resize(hot,extend,fix="center")  
  # get percentiles for scores
  mcols(hot)$perc=ecdf(hot$score)(hot$score)
  #' set chr
  bs.seqinfo = SeqinfoForUCSCGenome(bs)
  bs.seqinfo.commonchrs = keepSeqlevels(bs.seqinfo,seqlevels(hot) )
  seqinfo(hot)=bs.seqinfo.commonchrs
  #' remove short ones
  hot=hot[width(GenomicRanges::trim(hot)) == extend,]
  #' add percentiles
  mcols(hot)$perc=ecdf(hot$score)(hot$score)
  
  hot
}
hot=processHOTsummits_ucsc(
  "/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/HOThg19_filterNearBy2kb_Uniform_nopol2a_summits.bed",
  extend=2000, "hg19")
wind = add.hot.perc.intervals(hot)$hot

tmp = readRDS("/data/akalin/kwreczy/projects/HOT_DATA/rds/meth.hot.rds")
wind.meth = tmp[[2]]
cpgi.nonhot.meth = tmp[[3]]

wind.meth.hot = wind.meth[ which(wind$x=="0.99-1"), ]

wind.meth.hot.stats.col=data.frame(sample="0.99-1",
                                   sd=apply(wind.meth.hot, 2, sd, na.rm = TRUE),
                                   var=apply(wind.meth.hot, 2, var, na.rm = TRUE),
                                   ave=colMeans( wind.meth.hot, na.rm = TRUE ),
                                   med = apply(wind.meth.hot, 2, median, na.rm = TRUE),
                                   iqr=apply(wind.meth.hot, 2, IQR, na.rm = TRUE),
                                   stringsAsFactors = FALSE)
cpgi.nonhot.met.stats.col=data.frame(sample="non-HOT CpGi",
                                     sd=apply(cpgi.nonhot.meth, 2, sd, na.rm = TRUE),
                                     var=apply(cpgi.nonhot.meth, 2, var, na.rm = TRUE),
                                     ave=colMeans( cpgi.nonhot.meth , na.rm = TRUE),
                                     med = apply(cpgi.nonhot.meth, 2, median, na.rm = TRUE),
                                     iqr=apply(cpgi.nonhot.meth, 2, IQR, na.rm = TRUE),
                                     stringsAsFactors = FALSE)
stats.col = rbind(wind.meth.hot.stats.col, cpgi.nonhot.met.stats.col)
print(colnames(stats.col))
# for ggplot labels
colnames(stats.col) = c("sample", "Standard deviation" ,   "Variance"  ,"Average" ,   "Median"  ,  "Interquartile range")


library(cowplot)
dfmelt<-stats.col[,c("sample", "Median")]
colnames(dfmelt) <- c("sample", "value")
#dfmelt = cbind(dfmelt, cellline="Cell line") # add this if you want cell line in grey rectanfular above boxplots
#dfmelt$celline = as.character(dfmelt$cellline)# add this if you want cell line in grey rectanfular above boxplots
m_median<-ggplot(dfmelt, aes(x=sample, y=value))+
  #stat_boxplot(geom ='errorbar') + 
  geom_boxplot(colour="black", fill="steelblue", 
               outlier.size = 0.5)+
  #geom_boxplot(colour="black", fill="steelblue")+
  #facet_wrap(~ celline)+ # add this if you want cell line in grey rectanfular above boxplots
  labs(y="Median of methylation\nper cell line\n", x="") +
  scale_y_continuous(limits = c(0, 1))+
  my_theme 

dfmelt<-stats.col[,c("sample", "Interquartile range")]
colnames(dfmelt) <- c("sample", "value")
#dfmelt = cbind(dfmelt, cellline="Cell line")
#dfmelt$celline = as.character(dfmelt$cellline)
m_iqr<-ggplot(dfmelt, aes(x=sample, y=value))+
  #stat_boxplot(geom ='errorbar') + 
  geom_boxplot(colour="black", fill="steelblue", 
               outlier.size = 0.5)+
  #facet_wrap(~ celline)+
  #geom_boxplot(colour="black", fill="steelblue")+
  labs(y="Interquartile range of methylation\nper cell line\n", x="") +
  scale_y_continuous(limits = c(0, 1))+
  my_theme 

#print(tabl[6,]) #ESC.H9
i=6
df1 = data.frame(sample = c(as.character(wind$x), rep("non-HOT CpGi",nrow(cpgi.nonhot.meth))),
                 values = c(wind.meth[,i],cpgi.nonhot.meth[,i]),
                 #celline=tabl[,2][i],
                 celline="H9 cell line",
                 stringsAsFactors = F)
m_h9=ggplot(df1, aes(x = sample, y = values)) +
  #stat_boxplot(geom ='errorbar') + 
  geom_boxplot(colour="black", fill="steelblue", 
               outlier.size = 0.5)+
  #geom_boxplot(fill="steelblue") +
  facet_wrap(~ celline)+
  xlab("") + 
  ylab("Average methylation\n") +
  ggtitle("") +
  my_theme 


myplot1 = dripseq_GSE68948_GSE57353 
myplot2 = dripseq_GSE70189.1 
myplot3 = dripseq_input_ce10 
myplot4 = yeast_dat_m
myplot5 = g4chipseq_m

myplot6 = m_h9
myplot7 = m_median
myplot8 = m_iqr



white.space <- rectGrob(gp=gpar(fill="white"))
white.space1<- rectGrob(gp=gpar(fill="white"))
white.space2<- rectGrob(gp=gpar(fill="white"))

myplot1 <- arrangeGrob(myplot1, top = textGrob("A", x = unit(0, "npc")
                                                   , y   = unit(1, "npc"), just=c("left","top"),
                                                   gp=gpar(col="black", fontsize=12, fontfamily="Helvetica",fontface = "bold")))
myplot2 <- arrangeGrob(myplot2, top = textGrob("", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12, fontfamily="Helvetica",fontface = "bold")))
myplot3 <- arrangeGrob(myplot3, top = textGrob("B", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12, fontfamily="Helvetica",fontface = "bold")))
myplot4 <- arrangeGrob(myplot4, top = textGrob("C", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12, fontfamily="Helvetica",fontface = "bold")))
myplot5 <- arrangeGrob(myplot5, top = textGrob("D", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12, fontfamily="Helvetica",fontface = "bold")))
myplot6 <- arrangeGrob(myplot6, top = textGrob("E", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12, fontfamily="Helvetica",fontface = "bold")))
myplot7 <- arrangeGrob(myplot7, top = textGrob("F", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12, fontfamily="Helvetica",fontface = "bold")))
myplot8 <- arrangeGrob(myplot8, top = textGrob("", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12, fontfamily="Helvetica",fontface = "bold")))
myplotX <- arrangeGrob(white.space, top = textGrob("", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12, fontfamily="Helvetica",fontface = "bold")))


pdf("~/Figure4.pdf", width = 13, height =  7)
library(gridExtra)

grid.arrange(myplot1, myplot2, myplot3, myplot4,
             myplot5 ,myplot6,myplot7,myplot8,
             white.space1, white.space2,
             ncol=1,
             layout_matrix = rbind(c(1,1,1,1,2,2,3,3,4),
                                   c(5,5,5,9,6,6,7,8,10)))

dev.off()




