#'---
#'title: "Machine-learning models for HOT regions"
#'author: "Katarzyna Wreczycka"
#'output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    fig_width: 10
#'    fig_height: 10
#'---

setwd("/home/kwreczy/projects/HOT-regions/Generic_ML")
data_path="/data/akalin/Projects/AAkalin_HOTRegions/Data/"
output_path = "/data/akalin/kasia/projects/HOT_DATA/generic_ML_for_kmer_hot/"

library(GenomicRanges)
library(glmnet)
library(cowplot)

#######
##### FUNCTIONS
######

source("./hot_functions.R")
source("./MLreport_cpgi_functions.R")

#######
##### MAIN
######

#' # Goal
#' We built a predictive model of "hotness" of genomic regions using a penalized multivariate regression
#' method called Elastic Net for cross-species comparison of feature conservation in the sequence level.
#' Model variables includes observed GC content, observed
#' to expected GC ratio, GC skew and 2-,3-,4 mers. Elastic net method was performed not only within species,
#' but also across species. 
#'
#' O/E ratio is very predictive for mm9 and hg19, that's why in the second 
#' part of this script we additionally took it into account when training.
#'
#' We choose percentiles to decide which regions are hot:
#' dm3, ce10, mm9 = [0.99, 1],
#' hg19 = [0.995, 1]
#' "Cold" regions are with percentiles [0, 0.85].
#' 

#' # Create training and test data, run Elastic net for three different training sets

source("./generic_ML_for_kmer_hot.regions.R")
# It takes a while..
# paths.dm3 = createTrainingData.runElasticNet("dm3")
# paths.ce10 = createTrainingData.runElasticNet("ce10")
# paths.hg19 = createTrainingData.runElasticNet("hg19")
# paths.mm9 = createTrainingData.runElasticNet("mm9")
# paths.hg19 = createTrainingData.runElasticNet("hg19", cpgi_sampling=FALSE)
# paths.mm9 = createTrainingData.runElasticNet("mm9", cpgi_sampling=FALSE)

paths.dm3 = list(
  modelpath="/data/akalin/kwreczy/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelsdm3_nogcsk.rds",
  testpath= "/data/akalin/kwreczy/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.trainingdm3_nogcsk.rds")
paths.ce10 = list(
  modelpath ="/data/akalin/kwreczy/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelsce10_nogcsk.rds",
  testpath="/data/akalin/kwreczy/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.trainingce10_nogcsk.rds"
)
paths.hg19 = list(
  modelpath="/data/akalin/kwreczy/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelshg190.995_filterNearBy2kb_Uniform_nopol2a_nogcsk.rds",
  testpath="/data/akalin/kwreczy/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.traininghg19_filterNearBy2kb_Uniform_nopol2a_nogcsk.rds")
paths.mm9 = list(
  modelpath="/data/akalin/kwreczy/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelsmm9_nopol2a_nogcsk.rds",
  testpath="/data/akalin/kwreczy/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.trainingmm9_nopol2a_nogcsk.rds")

mods.path = list(hg19=paths.hg19$modelpath,
                 mm9=paths.mm9$modelpath,
                 dm3=paths.dm3$modelpath,
                 ce10=paths.ce10$modelpath)
test.path = list(hg19=paths.hg19$testpath,
                 mm9=paths.mm9$testpath,
                 dm3=paths.dm3$testpath,
                 ce10=paths.ce10$testpath)


org = c("hg19","mm9","dm3","ce10")
# All combination of 2 organisms
comb = combn(org, 2)
comb = data.frame(c(org, comb[1,],comb[2,]), 
                  c(org, comb[2,],comb[1,]),
                  stringsAsFactors = FALSE)
colnames(comb) = c("test","model")


#' # Models with coefficients - gcsk, obs, oe, GC skew, 2,3,4-mers and with decrease number of human HOT regions (filtering area 2kb)
#' ## Plot means of coefficients as importance

plot.aveimportance(org, mods.path)
plot.predscores.boxplots(comb, mods.path, test.path)

#' ## Heatmap showing AUC of different combinations of models and test data

#' Rows correspond to test data, columns to models
#ma = calc.AUC.matrix(org, mods.path, test.path)
#saveRDS(ma, paste0(output_path,"AUC.matrix_uniform_nopol2a.rds")) # filtering 1kb in hg19, no cpgi improvement in hg19 and mm9
ma= readRDS(paste0(output_path,"AUC.matrix_uniform_nopol2a.rds"))

plot.hist = function(ma, xlab, ylab, col){
  require(reshape2)
  require(ggplot2)
  require(RColorBrewer)
  
  melted_cormat <- melt(ma)
  g = ggplot(data = melted_cormat, aes(x=X2, y=X1, fill=value)) + 
    #theme_bw()+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=15))+
    geom_tile(aes(fill = value)) +
    geom_text(aes(fill = melted_cormat$value, label = round(melted_cormat$value, 2)), size=3) +
    scale_fill_gradientn(colours=col, name="AUC")+ 
    xlab(xlab) +
    ylab(ylab)
  #g+theme(axis.text=element_text(size=21),
  #       axis.title=element_text(size=21))
  g
}

pdf("/home/kwreczy/projects/HOT-regions/HOT-paper/HOT-paper-figures/AUCmatrix.pdf") #it's old path, I didn't save it anywhere yet
plot.hist(ma,
          ylab="Models for different species\n",
          xlab="\nSpecies specific test data",
          col=brewer.pal(9,"Purples")[5:8]) #without the darkest colour, with it text inside cells is not readable
dev.off()

#pdf("/home/kwreczy/projects/HOT-regions/HOT-paper/HOT-paper-figures/barplots_ml.pdf") #it's old path
plot.betas.barplots(org, mods.path)
#dev.off()

#' # Models are built with CpGi overlap based sampling for hg19 and mm9 and with decrease number of human HOT regions (filtering area 2kb). Models with coefficients - gcsk, obs, oe, GC skew, 2,3,4-mers
#' 

paths.hg19 = list(
  modelpath="/data/akalin/kasia/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelshg190.995_filterNearBy2kb_cpgi_Uniform_nopol2a_nogcsk.rds",
  testpath="/data/akalin/kasia/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.traininghg19_filterNearBy2kb_Uniform_nopol2a_nogcsk.rds")
paths.mm9 = list(
  modelpath="/data/akalin/kasia/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelsmm9_cpgi_nopol2a_nogcsk.rds",
  testpath="/data/akalin/kasia/projects/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.trainingmm9_nopol2a_nogcsk.rds")

mods.path = list(hg19=paths.hg19$modelpath,
                 mm9=paths.mm9$modelpath,
                 dm3=paths.dm3$modelpath,
                 ce10=paths.ce10$modelpath)
test.path = list(hg19=paths.hg19$testpath,
                 mm9=paths.mm9$testpath,
                 dm3=paths.dm3$testpath,
                 ce10=paths.ce10$testpath)
#' ## Plot means of coefficients as importance for calc. models
org1 = c("hg19", "mm9")
plot.aveimportance(org1, mods.path)
plot.predscores.boxplots(comb, mods.path, test.path)

#' ## Heatmap showing AUC of different combinations of models and test data
#' Rows correspond to test data, columns to models
#ma = calc.AUC.matrix(org, mods.path, test.path)
#saveRDS(ma, paste0(output_path,"AUC.matrix_cpgi_uniform_nopol2a.rds")) # filtering 2kb in hg19, cpgi improvement in hg19 and mm9
ma = readRDS(paste0(output_path,"AUC.matrix_cpgi_uniform_nopol2a.rds")) # filtering 1kb in hg19, no cpgi improvement in hg19 and mm9

#pdf("/home/kwreczy/projects/HOT-regions/HOT-paper/HOT-paper-figures/AUCmatrix_cpgi.pdf", width = 8, height = 6)
plot.hist(ma,
          ylab="Models for different species\n",
          xlab="\nSpecies specific test data",
          col=brewer.pal(9,"Purples")[4:8]) #wihtout the darkest colour, with it values inside cells are not readable
#dev.off()


#' ## Barplots showing normalized per species importance scores of features to 0-100 scale and average across species
#pdf("/home/kwreczy/projects/HOT-regions/HOT-paper/HOT-paper-figures/barplots_ml_cpgi.pdf")
plot.betas.barplots(org, mods.path)
#dev.off()

