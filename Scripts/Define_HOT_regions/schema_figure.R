
library(GenomicRanges)
# this object corresponds to TFBS from ENC TF Binding Super-track
tfs=readRDS("~/HOT/scripts/define_HOT_regions/data/tf_grl.rds") 
# this object corresponds to TFBS from Uniform TFBS Track
#tfs=readRDS("/data/akalin/kasia/projects/HOT_DATA/rds/human_tfbs_uniform.rds") # take this object, propably first use unlist() on it
r=coverage(tfs)

# get all peaks in one GRanges object
peaks=unlist(tfs)
end(peaks)=start(peaks) + peaks$peak
start(peaks)=end(peaks)

pcov500=runsum(pcov, k = 501,endrule="constant")

# get peak summit coverage 501 bp sliding window density
source("./peaker_adapted.R")
peakAnchorsW05kFl=getPeakRle(pcov500,filterNearBy=2000)


#' plots peaks, density of peaks and HOT regions anchors
#' 
#' @param peaks Chip-seq peaks in GRanges format
#' @param dens density of peak summits in RleList format
#' @param summits anchors for hot Regions in GRanges format
#' @param range a GRanges object with length 1 that defines where to center the plot
peakPlot4<-function(peaks,dens,summits,hotregions, 
                    range=NULL,chr=NULL,start=NULL,end=NULL,extend=500000,
                    col.peaks="black",col.dens="black", col.anchors="black", col.hot="black",...){
  require(Gviz)
  if(! is.null(range)){
    chr=as.character(seqnames(range[1,]))
    start=(start(range[1,]))
    end=(end(range[1,]))
  }
  wind=GRanges(chr,IRanges(start=start-extend,end=end+extend))
  
  peaks=AnnotationTrack(peaks[peaks %over% wind ,],name="peaks", col=col.peaks)
  
  dens=DataTrack(range=as(dens[chr],"GRanges"),chromosome=chr,type="histogram",
                 name="summitDensity",col.histogram=col.dens,
                 fill.histogram=col.dens)
  
  anchors = summits[summits %over% wind ,]
  summits=AnnotationTrack(anchors,name="Anchors", col=col.anchors)
  
  hotre=AnnotationTrack(hotregions, name="HOT regions", fill=col.hot)
  
  plotTracks(list( peaks, dens, summits, hotre,
                   GenomeAxisTrack()),
             from=start-extend, to=end+extend,...)
}


svg("~/HOT/scripts/define_HOT_regions/schema7.svg")
peakPlot4(peaks=unlist(tfs),dens=(pcov500),summits=peakAnchorsW05kFl,
          hotregions=hotreg[hotreg$perc >= 0.99],
          range=peakAnchorsW05kFl[(peakAnchorsW05kFl$scores>50),][4,],
          extend=10000,sizes=c(1.2,0.7,0.3,0.3, 0.3),
          background.panel = "white", background.title = "darkblue",
          col.peaks="black", col.dens="grey", col.anchors="black", col.hot="#FFB85F")
dev.off()


############ add genes annotation to figure
library(Gviz)
options(ucscChromosomeNames=FALSE)
library(biomaRt)
mart = useMart("ensembl",dataset="rnorvegicus_gene_ensembl")
myChr = 2
myStart = 230947473
myEnd = 230958234
# First create track with gene model:
biomTrack = BiomartGeneRegionTrack(genome="rn5", biomart=mart, chromosome=myChr,
                                   start=myStart, end=myEnd,showId=F, geneSymbols=F,
                                   rotate.title=TRUE, col.line=NULL, col="orange",
                                   fill="orange",filters=list(biotype="protein_coding"),
                                   collapseTranscripts=FALSE)
hotreg = resize( peakAnchorsW05kFl, 2000,fix="center")
mcols(hotreg)$perc=ecdf(hotreg$scores)(hotreg$scores)
peakPlot4(peaks=unlist(tfs),dens=(pcov500),summits=peakAnchorsW05kFl,
          hotregions=hotreg[hotreg$perc >= 0.99],
          range=peakAnchorsW05kFl[(peakAnchorsW05kFl$scores>50),][4,],
          extend=10000,sizes=c(1.2,0.7,0.3,0.3, 0.3),
          background.panel = "white", background.title = "darkblue")


