
.libPaths(c("/home/kwreczy/Rlibs/3.4/"))
data_path="/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/"

#######
##### FUNCTIONS
######


require(GenomicRanges)
require(Biostrings)
require(genomation)


#' get OE ratio and GC content for a given set of DNAstrings
GCinfo<-function(bs,wins, str.set)
{
  di.mat=dinucleotideFrequency(  str.set )
  a.mat =alphabetFrequency( str.set ,baseOnly=TRUE )
  gcsk=(a.mat[,3]-a.mat[,2])/(a.mat[,3] + a.mat[,2])
  exp=(a.mat[,2]*a.mat[,3])/width(str.set )
  obs=di.mat[,"CG"]
  OE=obs/exp
  
  cbind(gcsk=gcsk,
        obs=obs/width(str.set),
        oe=OE,
        gc=100*(a.mat[,2]+a.mat[,3])/rowSums(a.mat[,1:4]) #gc content
  )
  
}

#' process HOT regions and control GRanges objects
processHOTsummits<-function(fname,extend=1000,bs,add.chr=FALSE){
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
  seqinfo(hot)=SeqinfoForBSGenome(bs)
  
  #' remove short ones
  hot=hot[width(GenomicRanges::trim(hot)) == extend,]
  
  hot
}

skew<-function(str){
  x=letterFrequencyInSlidingView(str, view.width=1, letters=c("C","G"), OR="|", as.prob=FALSE)
  #sk=(x[,2]-x[,1])/(x[,2]+x[,1])
  #max(sk)-min(sk)
  G=-1*x[,1]+x[,2]
  G=cumsum(G)
  cbind(diff.skew=max(G)-min(G), max.skew=max(abs(G)) )
}
skew.paral<-function(i, str){
  str = str[[i]]
  x=letterFrequencyInSlidingView(str, view.width=1, letters=c("C","G"), OR="|", as.prob=FALSE)
  #sk=(x[,2]-x[,1])/(x[,2]+x[,1])
  #max(sk)-min(sk)
  G=-1*x[,1]+x[,2]
  G=cumsum(G)
  cbind(diff.skew=max(G)-min(G), max.skew=max(abs(G)) )
}

#######
##### MAIN
######


# read HOT regions from human
# process human regions
library("BSgenome.Hsapiens.UCSC.hg19")
# path to bed file with human HOT regions
fname=paste0(data_path, "HOThg19_filterNearBy2kb_Uniform_nopol2a_summits.bed")
# read HOT regions
hhot.hg19=processHOTsummits(fname=fname,
                            extend=1000,
                            Hsapiens,
                            add.chr=FALSE)  

# get percentiles for scores
mcols(hhot.hg19)$perc=ecdf(hhot.hg19$score)(hhot.hg19$score)

seqinfo(hhot.hg19)=SeqinfoForBSGenome(Hsapiens)
chrom=c(paste0("chr",1:22, sep=""), "chrX", "chrY")
keepSeqlevels(hhot.hg19, chrom)
trim(hhot.hg19)

#hhot.hg19=hhot.hg19[width(trim(hhot.hg19)) == 1000,]

# remove regions with less than 5 peaks
#hhot.hg19=hhot.hg19[hhot.hg19$score >7 ,]

#' get seqs
seqs=getSeq(Hsapiens,hhot.hg19)

#' get GC skew 
tmp = mclapply(1:length(seqs), skew.paral, str=seqs, mc.cores=50)
skews = do.call(rbind, tmp)

# get GC info
gcsk=GCinfo(Hsapiens,hhot.hg19, seqs)


.jets<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


pdf("~/cpgfreq_gccontent_hg19.pdf", width = 7, height = 5)
x=gcsk[,3]
y=gcsk[,4]
smoothScatter(x, y,
              colramp = .jets,
              main="",
              ylab="GC content", xlab="CpG frequency")
points(
  x[which(hhot.hg19$perc>=.99)],
  y[which(hhot.hg19$perc>=.99)],
  pch=20,
  col=rgb(0,0,0,alpha=0.1)
)
dev.off()



pdf("~/gcskew_gccontent_hg19.pdf", width = 7, height = 5)

x=skews[,1]
y=gcsk[,3]
smoothScatter(x, y,
              colramp = .jets,
              main="",
              xlab="GC skew", ylab="GC content")
points(
  x[which(hhot.hg19$perc>=.99)],
  y[which(hhot.hg19$perc>=.99)],
  pch=20,
  col=rgb(0,0,0,alpha=0.1)
)
dev.off()



# read HOT regions from mouse
# process mouse regions
library("BSgenome.Mmusculus.UCSC.mm9")
fname=paste0(data_path, "/HOTmm9_nopol2a_summits.bed")
hhot.hg19=processHOTsummits(fname=fname,
                            extend=1000,
                            Mmusculus,
                            add.chr=FALSE)  

# get percentiles for scores
mcols(hhot.hg19)$perc=ecdf(hhot.hg19$score)(hhot.hg19$score)

seqinfo(hhot.hg19)=SeqinfoForBSGenome(Mmusculus)
chrom=c(paste0("chr",1:22, sep=""), "chrX", "chrY")
keepSeqlevels(hhot.hg19, chrom)
hhot.hg19 <- trim(hhot.hg19)

hhot.hg19=hhot.hg19[width(trim(hhot.hg19)) == 1000,]

# remove regions with less than 5 peaks
#hhot.hg19=hhot.hg19[hhot.hg19$score > 5 ,]

#' get seqs
seqs=getSeq(Mmusculus,hhot.hg19)

#' get GC skew 
tmp = mclapply(1:length(seqs), skew.paral, str=seqs, mc.cores=50)
skews = do.call(rbind, tmp)

# get GC info
gcsk=GCinfo(Mmusculus,hhot.hg19, seqs)


.jets<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


pdf("~/cpgfreq_gccontent_mm9.pdf", width = 7, height = 5)
x=gcsk[,3]
y=gcsk[,4]
smoothScatter(x, y,
              colramp = .jets,
              main="",
              xlab="CpG frequency", ylab="GC content")
points(
  x[which(hhot.hg19$perc>=.99)],
  y[which(hhot.hg19$perc>=.99)],
  pch=20,
  col=rgb(0,0,0,alpha=0.1)
)
dev.off()


pdf("~/gcskew_gccontent_mm9.pdf", width = 7, height = 5)

x=skews[,1]
y=gcsk[,3]
smoothScatter(x, y,
              colramp = .jets,
              main="",
              xlab="GC skew", ylab="GC content")
points(
  x[which(hhot.hg19$perc>=.99)],
  y[which(hhot.hg19$perc>=.99)],
  pch=20,
  col=rgb(0,0,0,alpha=0.1)
)
dev.off()

