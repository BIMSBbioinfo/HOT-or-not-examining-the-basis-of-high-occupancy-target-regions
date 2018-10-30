#' # Check expression profiles of HOT region genes
#' https://www.dropbox.com/sh/hl2shd48urp53ip/AABVJLT2NxXPeM3IIYVkBFHYa?dl=0

#' read genes
hot=read.table("/home/kwreczy/projects/HOT-regions/expression_on_HOT/20160914-public-3.0.0-FwVjFU-mm9-all-gene.txt",header=F,sep="\t")
#' read expression profiles
f5Cellt=read.table("/home/kwreczy/projects/HOT-regions/expression_on_HOT/E-MTAB-3578-query-results.tsv",sep="\t",header=T)

#' #' read genes
#' hot=read.table("data/top99p2gene.byGREAT-mm9-all-gene.txt",header=F,sep="\t")
#' #' read expression profiles
#' f5Cellt=read.table("data/E-MTAB-3578-query-results.tsv",sep="\t",header=T)

#' change the NA expression profiles to 0
f5CelltnoNA=f5Cellt
f5CelltnoNA[is.na(f5CelltnoNA)]=0

#' ## get percentiles for each profile


#' change expression values to percentiles.
#' for the data with NA values changed to 0.
#' 
f5CelltnoNAp=f5CelltnoNA
for(i in 3:ncol(f5CelltnoNA)){
  f5CelltnoNAp[,i]=100*ecdf(f5CelltnoNAp[,i])(f5CelltnoNAp[,i])
}


#' change expression values to percentiles.
#' for the data with NA values
f5Celltp=f5Cellt
for(i in 3:ncol(f5Cellt)){
  f5Celltp[,i]=100*ecdf(f5Celltp[,i])(f5Celltp[,i])
}


#' ## sub set by hot genes and do boxplots
#' for the data with NA values
boxplot(f5Celltp[f5Celltp$Gene.Name %in% hot$V1, 3:ncol(f5Celltp)])

par(mar=c(4, 17, 3, 4))
boxplot(f5CelltnoNAp[f5CelltnoNAp$Gene.Name %in% hot$V1, 3:ncol(f5CelltnoNAp)],
        xlab="Expression percentile",
        #border=par("bg"),
        main="HOT-region genes expression per cell type",
        col="steelblue", 
        horizontal = TRUE,
        yaxt="n",
        ylim=c(0,100))

labs <-  colnames(f5Celltp)[3:ncol(f5CelltnoNAp)]
labs=gsub("\\."," ",labs)
axis(side=2,at=1:35,labels=labs,las=2,cex.axis=0.7)





