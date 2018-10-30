options(scipen=999) #disable scientific notation
#.libPaths(c("/home/kwreczy/Rlibs/3.4/"))

################
############ IMPORTS
###############

library(GenomicRanges)
library(genomation)
library(tools)
library(rtracklayer)
library(parallel)

################
############ MAIN
###############

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
                  header=T, sep="\t",stringsAsFactors = FALSE)
tf.path = "/data/akalin/Projects/AAkalin_HOTRegions/Data/wgEncodeAwgTfbsUniform/narrowPeak/"
l = paste0(tf.path, paste0(meta$tableName, ".narrowPeak.gz"))
tf_list_of_gr <- mclapply(l, bednarropeak2gr, mc.cores=20)
tf_grl <- GRangesList(tf_list_of_gr)
names(tf_grl) <- meta$Factor

# remove pol2a
polymerase = c("POLR3G","POLR2A")
meta.pol2a.inx = which(meta$Factor %in% polymerase)
tf_grl.nopol2a = tf_grl[-meta.pol2a.inx]
META = meta[-meta.pol2a.inx,]
ENCODE_GRL = tf_grl.nopol2a

# combine peaks from the same TF but different cell lines/labs
TF_names_all = names(ENCODE_GRL)
TF_names = unique(names(ENCODE_GRL)) # length(TF_names_unique) 161

# How may cell lines per TF?

how.many.cl.per.TF = mclapply(1:length(TF_names), function(i){
  META.TFi = META[which(TF_names_all %in% TF_names[i] ), ]
  ENCODE_GRL.TFi = ENCODE_GRL[ which(TF_names_all %in% TF_names[i] ) ]
  length(unique(META.TFi$Cell.Line)) 
  #sort(table(META.TFi$Cell.Line)) 
}, mc.cores=20)
how.many.cl.per.TF = unlist(how.many.cl.per.TF)
names(how.many.cl.per.TF) = TF_names
sort(how.many.cl.per.TF, decreasing = TRUE)

# TF can have many the same cl because different 
# labs genertd chip-seq for the same tF and the same cell line
cl.per.TF.grl = mclapply(1:length(TF_names), function(i){
  
  # GRL per TF
  META.TFi = META[which(TF_names_all %in% TF_names[i] ), ]
  ENCODE_GRL.TFi = ENCODE_GRL[ which(TF_names_all %in% TF_names[i] ) ]
  
  # GRL per TF per CL
  cls.tf.tbl = sort(table(META.TFi$Cell.Line)) 
  cls.tf = names(cls.tf.tbl)
  
  # combine TF that are from the same cell line into one GRanges object
  cls.tf.combined = lapply(cls.tf, function(cl){
    unlist(ENCODE_GRL.TFi[ which(META.TFi$Cell.Line %in% cl) ])
  })
  names(cls.tf.combined) = cls.tf
  
  return(cls.tf.combined)
  
  
}, mc.cores=20)
names(cl.per.TF.grl) = TF_names

sapply(cl.per.TF.grl, length)



TFs.grl = mclapply(1:length(TF_names), function(i){
  
  # GRL per TF
  META.TFi = META[which(TF_names_all %in% TF_names[i] ), ]
  ENCODE_GRL.TFi = ENCODE_GRL[ which(TF_names_all %in% TF_names[i] ) ]
  
  # combine different cell lines, different labs together
  # into one GRanges per TF
  unlist(ENCODE_GRL.TFi)
  
}, mc.cores=20)
names(TFs.grl) = TF_names


TFs.grl.quality = mclapply(1:length(TF_names), function(i){
  
  # GRL per TF
  META.TFi = META[which(TF_names_all %in% TF_names[i] ), ]
  data.frame(caution=sum(META.TFi$quality=="caution"),
             good=sum(META.TFi$quality=="good"))
  
}, mc.cores=20)

TFs.quality = do.call("rbind", TFs.grl.quality)
rownames(TFs.quality) = TF_names
TF.perc.good.quality = TFs.quality$good / rowSums(TFs.quality) *100


source("~/projects/HOT-regions/basic_functions.R")
#BS.genomeSpecies=get.BS.genomeSpecies("hg19")
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
wind = keepSeqlevels(wind, unique(as.character(seqnames(wind))) )
wind.inter = split(wind, wind$x)


##################################################################################################
# Number of TFs per HOT region

calculate.nbr.Tfs.on.region = function(TFs.grl, region){
  
  HOT=region
  HOT$id = 1:length(HOT)
  
  TFs.grl.1 = GRangesList(TFs.grl)
  fo = findOverlaps(HOT, TFs.grl.1)
  
  ma.1 = matrix(0, length(HOT), length(TFs.grl.1))
  for(i in 1:length(HOT)){
    print(i)
    fo.hot.i = fo[ which( queryHits(fo) %in% i ) ]
    ma.1[i, subjectHits(fo.hot.i)] = 1
    
  }
  colnames(ma.1) = names(TFs.grl)
  return(ma.1)
}

ma.TFsonHOT = calculate.nbr.Tfs.on.region(TFs.grl, wind.inter[["0.99-1"]])


pdf("~/projects/HOT-regions/Review_NAR/TFsonHOT1.pdf", width = 30, height = 5)
# Reorder matrix
ord=order(colSums(ma.TFsonHOT ),
          decreasing = TRUE)
ma.TFsonHOT.orderedSums = ma.TFsonHOT[,ord]
# Create annotation heatmaps
annot_df1 <- data.frame(Perc.good.quality.antibody=TF.perc.good.quality[ord])
annot_df1_col = list(Perc.good.quality.antibody=colorRamp2(  c(0, 100),  c("red", "green")))

ha.top <- HeatmapAnnotation(annot_df1,
                            col=annot_df1_col,
                            barplot2 = anno_barplot(colSums(ma.TFsonHOT.orderedSums ),
                                                    gp = gpar(fill = c( "cornflowerblue")),
                                                    axis = TRUE),
                            height = unit(2, "cm"))
ha.right = rowAnnotation(b2 = row_anno_points(rowSums(ma.TFsonHOT.orderedSums ),
                                              axis = TRUE,
                                              width= unit(4, "cm")),
                         width= unit(4, "cm"))
# Add a legend
lgd = Legend(at = c("Yes", "No"),
             title = "Presence of TFs on HOT regions",
             type="points",
             legend_gp = gpar(col = c("red","blue"))
)
# Combine heatmaps
ht = Heatmap(ma.TFsonHOT.orderedSums,
             show_heatmap_legend = FALSE,
             #name = "Presence of TFs on HOT regions",
             row_title = "HOT regions",
             column_title = "Transcription factors",
             cluster_rows=FALSE,
             cluster_columns=FALSE,
             top_annotation = ha.top
) + ha.right
draw(ht, annotation_legend_list = list(lgd))

dev.off()


