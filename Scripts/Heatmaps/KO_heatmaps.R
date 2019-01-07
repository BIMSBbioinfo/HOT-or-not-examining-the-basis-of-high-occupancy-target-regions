

#### IMPORTS

library(GenomicRanges)
library(Rsamtools)
library(genomation)
library(RCurl)

#### FUNCTIONS

source("./Define_HOT_regions/extend_summits_and_samplewind.R")

# Calculates IP/control for given bam file paths
enrichmentMatrix<-function(signal,background,windows,type="bam",
                           extend.signal=0,extend.background=0,
                           lib.size.signal=NA, lib.size.background=NA,
                           pseudo.cnt=1,norm=TRUE,
                           bin.num=NULL, ...){
  message("getting the signal and background... ",date())
  if(!is.null(bin.num)){
    print(paste0("binning using ", bin.num, " bins"))
    signl=ScoreMatrixBin(signal,windows,extend=extend.signal,type=type,bin.num=bin.num,...)
    bckgr=ScoreMatrixBin(background,windows,extend=extend.background,type=type,bin.num=bin.num,...)
  }else{
    signl=ScoreMatrix(signal,windows,extend=extend.signal,type=type,...)
    bckgr=ScoreMatrix(background,windows,extend=extend.background,type=type,...)
    
  }
  if(norm){
    message("normalizing... ",date())
    command=paste0("samtools idxstats ",signal," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
    sig.mapped=as.numeric(try(system(command,intern=TRUE)))
    command2=paste0("samtools idxstats ",background," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
    bg.mapped=as.numeric(try(system(command2,intern=TRUE)))
    
    signl=1e6*signl/sig.mapped
    bckgr=1e6*bckgr/bg.mapped
    
  }else{
    signl=1e6*signl/lib.size.signal
    bckgr=1e6*bckgr/lib.size.background
  }
  log2((signl+pseudo.cnt)/(bckgr+pseudo.cnt))
}

enrichmentMatrix_paral = function(i, signals, backgrounds, windows, type, 
                                  extend.signal, extend.background,
                                  lib.size.signal, lib.size.background, norm, bin.num=NULL){
  enrichmentMatrix(signals[i],backgrounds[i],windows,type=type,
                   extend.signal=extend.signal[i], extend.background=extend.background[i],
                   lib.size.signal=lib.size.signal[i], lib.size.background = lib.size.background[i], 
                   norm=norm, bin.num=bin.num)
}  

enrichmentMatrix.no.ctrl<-function(signal,windows,type="bam",
                                   extend.signal=0,
                                   lib.size.signal=NA, 
                                   pseudo.cnt=1,norm=TRUE,
                                   bin.num=NULL, ...){
  message("getting the signal and background... ",date())
  if(!is.null(bin.num)){
    print(paste0("binning using ", bin.num, " bins"))
    signl=ScoreMatrixBin(signal,windows,extend=extend.signal,type=type,bin.num=bin.num,...)
  }else{
    signl=ScoreMatrix(signal,windows,extend=extend.signal,type=type,...)
    
  }
  if(norm){
    message("normalizing... ",date())
    command=paste0("samtools idxstats ",signal," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
    sig.mapped=as.numeric(try(system(command,intern=TRUE)))
    
    signl=1e6*signl/sig.mapped
    
  }else{
    signl=1e6*signl/lib.size.signal
  }
  signl
}
enrichmentMatrix.no.ctrl_paral = function(i, signals, windows, type, 
                                          extend.signal, 
                                          lib.size.signal,  norm, bin.num=NULL){
  enrichmentMatrix.no.ctrl(signals[i],windows,type=type,
                           extend.signal=extend.signal[i], 
                           lib.size.signal=lib.size.signal[i], 
                           norm=norm, bin.num=bin.num)
}  


.jets<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
.jets2<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                           "#FF7F00", "red", "#7F0000"))


#' Read murine HOT regions
hot=readBed("HOTmm9_nopol2a_summits.bed")
hhot <- resize(hot,2000,fix="center")
mcols(hhot)$perc=ecdf(hhot$score)(hhot$score)
hhot.mm9 = hhot[order(-hhot$score)]
hhot.mm9$scores = hhot.mm9$score
set.seed(111)
wind.sampled = sample.windows(hhot.mm9, k=3000)



#' Get information of KO samples about their positive or negative association with HOT region score,
#' so information about in which KO experiments we can observe enrichment of TFBS on HOT regions.
library(RCurl)
link="https://docs.google.com/spreadsheets/d/1C5AST4gjOo7mngjB95zsXLXwmDoP9YqWIVFASMEonaU/pub?gid=1861042470&single=true&output=csv"
myCsv <- getURL(link)
ishot <- read.csv(textConnection(myCsv), header=TRUE, stringsAsFactors =FALSE)
hot.samples.sign.geo = ishot[which(ishot$is.hot==TRUE),]$sign.geo

#' Get info about location of KO Chip-seq bam files on disk
#' 
library(RCurl)
link2="https://docs.google.com/spreadsheets/d/1C5AST4gjOo7mngjB95zsXLXwmDoP9YqWIVFASMEonaU/pub?gid=617642946&single=true&output=csv"
myCsv2 <- getURL(link2)
sheet2 <- read.csv(textConnection(myCsv2), header=TRUE, stringsAsFactors =FALSE)


kos.with.ctrl =sheet2[-which( is.na(sheet2$bckgr.bam) ), ]

# sml.with.control <- mclapply(1:length(kos.with.ctrl$sign.bam),
#                 enrichmentMatrix_paral,
#                 signals=kos.with.ctrl$sign.bam,
#                 backgrounds=kos.with.ctrl$bckgr.bam,
#                 windows=wind.sampled,
#                 type="bam",
#                 extend.signal=rep(150,length(kos.with.ctrl$sign.bam)),
#                 extend.background=rep(150,length(kos.with.ctrl$sign.bam)),
#                 lib.size.signal=NA,
#                 lib.size.background=NA,
#                 norm=TRUE,
#                 bin.num=50,
#                 mc.cores=20)
# sml.with.control <- as( sml.with.control, "ScoreMatrixList")
# names(sml.with.control) <- kos.with.ctrl$name
# saveRDS(sml.with.control, "sml_with_ctrl_mm9nopol2a_binmun50.rds")
sml.with.control = readRDS("sml_with_ctrl_mm9nopol2a_binmun50.rds")


#### KO - positive associ. HOT regions

print(df.wc [ which(df.wc$diff > 0), ]$name)
name.peak =  c( "MLL4_d7"     ,       "MLL4_d2"     ,       "E2AKO1h"    ,      
                "E2AKO6h"     ,       "Stat6"     ,         "SRC-2"             ,
                "RAP1ko2"      ,      "runx1"     ,         "MLL4_myocytes"     ,
                "rxra"          ,     "MLL4_d0"   ,         "IRF4post"          ,
                "Stat4" , "RAP1ko1" )
sml.peak = sml.with.control[   match(name.peak,
                                     names(sml.with.control)) ]
##
## line plot
##
# par(mar=c(2,2,2,4), oma=c(2,2,2,5))
# plotMeta(sml.peak, profile.names = names(sml.peak), xcoords=c(-1000, 1000))

##
## heatmaps
##
multiHeatMatrix(sml.peak, legend.name =rep("log2(IP/control)", 
                                           length(sml.peak)), winsorize = c(0.5, 99),
                grid=FALSE,cex.lab=0.7, xcoords=c(-1000, 1000),
                matrix.main=names(sml.peak), group=lvs,
                common.scale=T,col=.jets2(7),
                cex.axis=0.5,
                cex.main=0.5)

##
## boxplots
##


#' ## Boxplots for all windows (very big variables)

# sml.my <- mclapply(1:length(sheet.myhot$sign.bam),
#                    enrichmentMatrix_paral,
#                 signals=sheet.myhot$sign.bam,
#                 backgrounds=sheet.myhot$bckgr.bam,
#                 windows=hhot.mm9,
#                 type="bam",
#                 extend.signal=rep(150,length(sheet.myhot$sign.bam)),
#                 extend.background=rep(150,length(sheet.myhot$sign.bam)),
#                 lib.size.signal=NA,
#                 lib.size.background=NA,
#                 norm=TRUE,
#                 bin.num=50,
#                 mc.cores=length(sheet.myhot$sign.bam))
# saveRDS(sml.my, "KO_positive_valuesforboxplots_all_hot_windows_nopol2a_M.sml.my.rds")
sml.my <-readRDS("KO_positive_valuesforboxplots_all_hot_windows_nopol2a_M.sml.my.rds")

# sml1 =  sml.my
# result = do.call(cbind, sml1)
# result.means = rowMeans(result, na.rm = TRUE) 
# saveRDS(result.means, "KO_positive_valuesforboxplots_allwindowsresult.ave.result.means.rds")
result.means = readRDS("KO_positive_valuesforboxplots_allwindowsresult.ave.result.means.rds")


# scorematrix removes windows out of chrs, so lets remove some windows from hot variable
# source("~/projects/HOT-regions/basic_functions.R")
# tfname="HOTmm9_nopol2a_summits.bed"
# hot=processHOTsummits_ucsc(fname=tfname,
#                       extend=2000,
#                       "mm9",
#                       add.chr=FALSE)
# mcols(hot)$perc=ecdf(hot$score)(hot$score)
# 
# hot$x = "0"
# hot[hot$perc > 0.99 &  hot$perc <=1]$x  = "0.99-1"
# hot[hot$perc > 0.9 &  hot$perc <=0.99]$x = "0.9-0.99"
# hot[hot$perc > 0.8 &  hot$perc <=0.9]$x = "0.8-0.9"
# hot[hot$perc > 0.7 &  hot$perc <=0.8]$x = "0.7-0.8"
# hot[hot$perc > 0.6 &  hot$perc <=0.7]$x = "0.6-0.7"
# hot[hot$perc > 0.5 &  hot$perc <=0.6]$x = "0.5-0.6"
# hot[hot$perc > 0.0 &  hot$perc <=0.5]$x = "0.0-0.5"
# hot$x = as.factor(hot$x) 
# 
# result1 = data.frame(x=as.numeric(result.means), 
#                      perc=as.numeric(hot$perc),
#                      y=as.character(hot$x))
#saveRDS(result1, "KO_positive_valuesforboxplots_allwindowsresult.ave.rds")

normalize01 = function(x) (x-min(x))/(max(x)-min(x))
result1$x1 = normalize01(result1$x)

require(ggplot2)
library(cowplot)
p <- ggplot(result1, aes(x=y, y=x1)) + 
  geom_boxplot(colour="black", fill="steelblue", outlier.size = 0.5)+
  labs(title="",x="HOT score percentiles", y = "Average log2(IP/control)")+
  theme(
    axis.text=element_text(size=5),
    axis.text.x = element_text(angle = 45,hjust = 1, size=9),
    axis.text.y = element_text(hjust = 1, size=9),
    axis.title=element_text(size=10,face="plain"),
    axis.line = element_line(size=1, colour = "black"), 
    panel.grid.major = element_line(colour = "#E6E6E6"))
print(p)  


#### KO - negative associ. HOT regions
name.valley = c( "NFI"      ,        "Tbet"        ,     "KAP1"    ,         "MLL4_EBPbeta_Cre",
                 "G9aGLP"    ,       "RAP1ko1"      ,    "SA1"      ,        "KAP1_KAP1"       ,
                 "Prdm16rep1" , "Prdm16rep2")
sml.valley = sml.with.control[ match(name.valley,
                                     names(sml.with.control)) ]
names(sml.valley) =    name.valley

##
## line plot
##
par(mar=c(2,2,2,4), oma=c(2,2,2,5))
plotMeta(sml.valley, profile.names = names(sml.valley), xcoords=c(-1000, 1000))

##
## heatmaps
##
multiHeatMatrix(sml.valley, legend.name =rep("", 
                                             length(sml.valley)), winsorize = c(0.5, 99),
                grid=FALSE,cex.lab=0.7, xcoords=c(-1000, 1000),
                matrix.main=names(sml.valley), group=lvs,
                common.scale=T,col=.jets2(7),
                cex.axis=0.5,
                cex.main=0.5)






