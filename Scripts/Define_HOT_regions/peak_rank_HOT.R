

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
#saveRDS(tf_grl, "/data/akalin/kasia/projects/HOT_DATA/rds/hg19_tfbs_uniform.rds")
#tf_grl <- readRDS("/data/akalin/kasia/projects/HOT_DATA/rds/hg19_tfbs_uniform.rds")

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


calc.matrix.median.values = function(wind.inter, TFs.grl, scale=TRUE){
  
  require(scales)
  
  matrix_median_scores = matrix(0, length(wind.inter), length(TFs.grl))
  
  median_scores = mclapply(1:length(wind.inter), function(w){
    
    regions = wind.inter[[w]] # HOT or non-HOT regions
    lapply(1:length(TFs.grl), function(j){
      
      tf = TFs.grl[[j]]
      tf <- tf [order(tf$signalValue, decreasing = TRUE),]

      if(scale){
        tf$rank_scaled = rescale(1:length(tf), to = c(0, 1000))
      }else{
        tf$rank_scaled = 1:length(tf)
      }

      tfpeaks.on.regions = subsetByOverlaps(tf, regions) 
      #out=median(tfpeaks.on.regions$signalValue)
      out=median(tfpeaks.on.regions$rank_scaled)
      out
    })
  }, mc.cores=20)
  
  for(row in 1:length(median_scores)){
    for(col in 1:length(median_scores[[1]])){
      matrix_median_scores[row, col] = median_scores[[row]][[col]]
    }
  }
  rownames(matrix_median_scores) = names(wind.inter)
  colnames(matrix_median_scores) = names(TFs.grl)
  
  return(matrix_median_scores)
}

matrix_median_rankes = calc.matrix.median.values(wind.inter, TFs.grl, scale=FALSE)
matrix_median_rankes_scaled = calc.matrix.median.values(wind.inter, TFs.grl)


pdf("~/peak_matrix_median_ranks_scaledBycol_orderedBySum.pdf", width = 35)
#pdf("~/peak_matrix_median_rankes_scaled_col_orderbyHOTvals.pdf", width = 30)
library(circlize)
library("ComplexHeatmap")
require(scales)
ma_complexheatmap = apply(matrix_median_rankes, 
                                        2, 
                                        function(x) rescale(x, to = c(0, 100)))


#ord = order( ma_complexheatmap[which(rownames(ma_complexheatmap)=="0.99-1"), ],decreasing = FALSE)
#ma_complexheatmap = ma_complexheatmap [ , ord]
ord=order( colSums(ma_complexheatmap) ,decreasing = TRUE)
ma_complexheatmap = ma_complexheatmap [ ,ord]

annot_df <- data.frame(TF.perc.good.quality.antibody=TF.perc.good.quality,
                       Number.cellines = how.many.cl.per.TF)
annot_df = annot_df[ord,]


number.celline.val = c(    1:10,  14 , 70 )
number.celline.col = c(rep("gold", 2),
           rep("darkgoldenrod3", 2),
           rep("darkolivegreen2", 2),
           rep("darkturquoise", 2),
           rep("deepskyblue4", 2),
           "darkorchid3",
           "darkorchid4")
names(number.celline.col) = number.celline.val
col = list(TF.perc.good.quality.antibody=colorRamp2(  c(0, 100),  c("red", "green")),
           Number.cellines=number.celline.col)

library(matrixStats)
ha.top <- HeatmapAnnotation(annot_df, 
                            barplot2 = anno_barplot(colSums(ma_complexheatmap),
                                                    gp = gpar(fill = c( "cornflowerblue")),
                                                    axis = TRUE),
                            col = col,
                            na_col = "grey")

df.ma_complexheatmap = as.matrix.data.frame(ma_complexheatmap, stringsAsFactors=FALSE)
df.ma_complexheatmap[is.na(df.ma_complexheatmap)] <- 100 
ha.right = rowAnnotation(b2 = row_anno_boxplot(df.ma_complexheatmap),
                         width = unit(3, "cm"))

Heatmap(ma_complexheatmap, 
        name = "Scaled median of peak ranks",
        row_title = "Regions",
        column_title = "Transcription factors",
        top_annotation = ha.top,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        na_col="grey"
        ) + ha.right

dev.off()



