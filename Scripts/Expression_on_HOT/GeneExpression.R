
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
## Regions
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
#wind = keepSeqlevels(wind, unique(as.character(seqnames(wind))) )
wind.inter = split(wind, wind$x)
#HOT (>99th percentile), MILD (between 99th and 75th percentile), and COLD regions (below 75th percentile)
HOT=wind[ wind$perc>.99 ]
MILD=wind[ wind$perc<.99 & wind$perc>.75 ]
COLD=wind[ wind$perc<.75 ]

# https://www.tau.ac.il/~elieis/HKG/
housekeeping.genes = read.table("/data/akalin/Projects/AAkalin_HOTRegions/Results/Reviews_NAR/HK_genes.txt", stringsAsFactors = FALSE)
library("biomaRt")
ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")
values<- housekeeping.genes[,2]
tbl.convert.genes = getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), filters = "refseq_mrna", values = values, mart= ensembl)

tbl.convert.housekeeping = tbl.convert.genes[
  match( tbl.convert.genes$refseq_mrna, housekeeping.genes[,2]),
  ]
housekeeping.genes.ensembl = tbl.convert.housekeeping$ensembl_gene_id #4231 h. genes

# Load ensembl genes
ensembl.genes = readBed( # 205000
  "/data/akalin/Projects/AAkalin_HOTRegions/Results/Reviews_NAR/ensembl_hg19.bed.gz")
ensembl.genes$genename = gsub("T", "G", ensembl.genes$name)

# Find Genes associated with HOT regions
library(rGREAT)
job = submitGreatJob("/data/akalin/Projects/AAkalin_HOTRegions/Results/HOT_regions_final/hg19-HOTregions.bed", species = "hg19")
tb = getEnrichmentTables(job)
res = plotRegionGeneAssociationGraphs(job)
res.2kb = res[ which(abs(res$distTSS) < 2000), ] # 3147


# Find ENSEMBL ids of Genes associated with HOT regions

#convert.ensemblToGeneName = read.table("~/ensembl_hg19_ensemblToGeneName.bed.gz", col.names = c("ensembl.transcriptname", "genename"))
library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# query biomart
results1 <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
                 #filters = "ensembl_transcript_id", 
                 #values = convert.ensemblToGeneName$ensembl.transcriptname,
                 filters="external_gene_name",
                 values=res.2kb$gene,
                 mart = mart)
results1 = results1[which(results1$gene_biotype == "protein_coding"),] # 3393
res.2kb.1 = res.2kb[ res.2kb$gene %in% intersect(res.2kb$gene, results1$external_gene_name) ]
tmp = results1[ match(res.2kb.1$gene, results1$external_gene_name ), ]
res.2kb.1$external_gene_name = tmp$external_gene_name
res.2kb.1$ensembl_gene_id =  tmp$ensembl_gene_id
  

HOT.genes = res.2kb.1
HOT_genes_ensembl = HOT.genes$ensembl_gene_id

#############################################################################
## rpkm values (ensembl ids)
#############################################################################

cellines= read.table("/data/akalin/Projects/AAkalin_HOTRegions/Data/RNAseq_hg19_ENCODE/RoadmapEpigenomics/EG.name.txt",
                     header=FALSE, 
                     sep="\t", 
                     stringsAsFactors = FALSE)
rpkm.matrix = read.table("/data/akalin/Projects/AAkalin_HOTRegions/Data/RNAseq_hg19_ENCODE/RoadmapEpigenomics/57epigenomes.RPKM.pc",  
                         header=TRUE,
                         sep="\t", 
                         stringsAsFactors = FALSE)
# fix names of columns
colnames(rpkm.matrix) = colnames(rpkm.matrix)[2:length(colnames(rpkm.matrix))]
rpkm.matrix=rpkm.matrix[,-length(colnames(rpkm.matrix)),drop=FALSE]
#match( cellines[,1], colnames(rpkm.matrix) ) == match(colnames(rpkm.matrix), cellines[,1] )
colnames(rpkm.matrix) = cellines[,2]


expression.housekeeping = rpkm.matrix[ which(rownames(rpkm.matrix) %in% housekeeping.genes.ensembl), ] #3239
expression.hot = rpkm.matrix[ which(rownames(rpkm.matrix) %in% HOT_genes_ensembl), ] # 2758
expression.nonhot = rpkm.matrix[ -which(rownames(rpkm.matrix) %in% HOT_genes_ensembl), ] # 17037



#############################################################################
## Plot it!
#############################################################################

# Just to double-check

# stolen from genomation
.winsorize<-function(mat,rng){
  hi.th=quantile(mat,rng[2]/100,na.rm=TRUE)
  lo.th=quantile(mat,rng[1]/100,na.rm=TRUE)
  mat[mat>hi.th]=hi.th
  mat[mat<lo.th]=lo.th
  mat
}
expression.hot.wind =.winsorize(expression.hot, c(1,90))
rescale <- function(x) (x-min(x))/(max(x) - min(x)) * 100
expression.hot.scaled = apply(expression.hot.wind ,2 , function(x) rescale(x))


# Plot
pdf("~/expression.hot.percentiles.hg19.pdf", width = 20, height=10)
op <- par(mar=c(25,4,1,1))
boxplot(expression.hot.scaled ,
        main="HOT-region gene expression per cell type/tissue", 
        xlab="", 
        ylab="Expression percentile",
        col = "cornflowerblue",
        las=2)

par(op)
dev.off()


g.jean = function(mymatrix){
  require(matrixStats)
  apply(as.matrix(mymatrix), 1, function(x) mad(x, na.rm = TRUE)/median(x, na.rm = TRUE))
}
expression.housekeeping.rowmed = g.jean(expression.housekeeping) 
expression.hot.rowmed  = g.jean(expression.hot)
expression.nonhot.rowmed =  g.jean(expression.nonhot)

bxpl1 = data.frame(var="housekeeping genes", value = expression.housekeeping.rowmed, stringsAsFactors = FALSE)
bxpl2 = data.frame(var="HOT genes", value = expression.hot.rowmed, stringsAsFactors = FALSE)
bxpl3 = data.frame(var="non-HOT genes", value = expression.nonhot.rowmed, stringsAsFactors = FALSE)
bxpl = rbind(rbind(bxpl1, means.bxpl2),bxpl3)
bxpl$var <- factor(bxpl$var,
                       levels = c('HOT genes','housekeeping genes','non-HOT genes'),ordered = TRUE)


library(genomation)
library(Rsamtools)
library(ggplot2)
library("gridExtra")
library(reshape)
library(cowplot)
library(plyr)
library(reshape)

pdf("~/gene_expression_HOT.pdf")
require(ggplot2)
ggplot(data = bxpl, aes(x=var, y=value)) + 
  geom_boxplot(aes(fill=var), show.legend=FALSE)+
  labs(x="Regions", y="Gene expression [MAD/median]", colour="") +
  scale_fill_manual(values=c("#FF420E", "#A2C523", "#336B87"))+
  theme(
    axis.text=element_text(size=10),
    axis.text.x = element_text(angle=45,hjust = 1, size=15),
    axis.text.y = element_text(hjust = 1, size=15),
    strip.text.x = element_text(size=20),
    axis.title=element_text(size=25, face="plain"),
    axis.line = element_line(size=1, colour = "black"), 
    panel.grid.major = element_line(colour = "#E6E6E6"))
dev.off()






