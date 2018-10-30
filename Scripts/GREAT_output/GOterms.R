 
#' # GO Biological Processes associated with HOT regions


#' ## mm9

tsv.data <- read.delim("./mm9/shown-GOBiologicalProcess.tsv",
                        header=T, sep="\t", skip=1)

tsv.data$x = -log10(tsv.data$Binom.Raw.P.Value)
sortd = tsv.data[with(tsv.data, order(-x)), ]
rownames(sortd)
# to break line in some rownames
rownames(tsv.data)[9] <- "intrinsic apoptotic signaling pathway\nin response to DNA damage"                      
rownames(tsv.data)[10] <- "regulation of cyclin-dependent protein\nserine/threonine kinase activity"
rownames(tsv.data)[14] <- "intrinsic apoptotic signaling pathway in response\nto DNA damage by p53 class mediator"

pdf("./mm9/GEOprocessesbarplot_mm9.pdf")
parOrg=par(mar=c(5, 4, 4, 2) + 0.1)
par(mar=c(4, 15, 1, 0) + 0.1)
barplot(rev(tsv.data$x), horiz=TRUE, col="steelblue", 
        cex.names=0.70,
        names.arg=rev(rownames(tsv.data)),las=1,
        cex.axis=1,
        xlab="-log10(Binomial p value)", main="GO Biological Process")
dev.off()

#' ## hg19

tsv.data <- read.delim("./hg19/shown-GOBiologicalProcess.tsv",
                       header=T, sep="\t", skip=1)

tsv.data$x = -log10(tsv.data$Binom.Raw.P.Value)
sortd = tsv.data[with(tsv.data, order(-x)), ]
rownames(sortd)

pdf("./hg19/GEOprocessesbarplot_hg19.pdf")
parOrg=par(mar=c(5, 4, 4, 2) + 0.1)
par(mar=c(4, 18, 1, 0.5) + 0.1)
barplot(rev(tsv.data$x), horiz=TRUE, col="steelblue", 
        cex.names=0.70,
        names.arg=rev(rownames(tsv.data)),las=1,
        cex.axis=1,
        xlab="-log10(Binomial p value)", main="GO Biological Process")
dev.off()





