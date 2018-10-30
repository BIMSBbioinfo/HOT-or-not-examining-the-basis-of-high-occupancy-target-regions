#'---
#'title: "PCA of HOT and COLD regions"
#'author: "Katarzyna Wreczycka"
#'output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    fig_width: 10
#'    fig_height: 10
#'---
#'


paths.dm3 = list(
    modelpath="/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelsdm3_nogcsk.rds",
    testpath= "/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.trainingdm3_nogcsk.rds")
paths.ce10 = list(
    modelpath ="/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelsce10_nogcsk.rds",
    testpath="/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.trainingce10_nogcsk.rds"
)
paths.hg19 = list(
    modelpath="/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelshg190.995_filterNearBy2kb_Uniform_nopol2a_nogcsk.rds",
    testpath="/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.traininghg19_filterNearBy2kb_Uniform_nopol2a_nogcsk.rds")
paths.mm9 = list(
    modelpath="/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.sampled.ML.modelsmm9_nopol2a_nogcsk.rds",
    testpath="/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/HOT.kmer.trainingmm9_nopol2a_nogcsk.rds")

mods.path = list(hg19=paths.hg19$modelpath,
                 mm9=paths.mm9$modelpath,
                 dm3=paths.dm3$modelpath,
                 ce10=paths.ce10$modelpath)
test.path = list(hg19=paths.hg19$testpath,
                 mm9=paths.mm9$testpath,
                 dm3=paths.dm3$testpath,
                 ce10=paths.ce10$testpath)


hg19 = readRDS(test.path$hg19)
dm3 = readRDS(test.path$dm3)
ce10 = readRDS(test.path$ce10)
mm9 = readRDS(test.path$mm9)
orgs = list(hg19, mm9, dm3, ce10)
names(orgs) = c("hg19", "mm9", "dm3","ce10")
percs = lapply(orgs, function(x) x[,1])
names(percs) = c("hg19", "mm9","dm3", "ce10")

#' HOT regions are >=.99 percentile and COLD regions <.50 percentile

data = lapply(names(orgs), function(name){
  print(name)
  perc = percs[[name]]
  which.hot = which(perc>=.99)
  print(length(which.hot))
  which.cold = which(perc<.5)
  which.cold.sampled =  sample(which.cold, length(which.hot))
  which.perc = c(which.hot, which.cold.sampled )
  list(
    len=length(which.perc),
    data=orgs[[name]][which.perc,]
  )
})
len.perc.orgs = sapply(data, function(x) x$len)
data = do.call(rbind, lapply(data, function(x) x$data))

#' # PCA (with center = TRUE and scale = TRUE)
# https://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
data.pca <- prcomp(data[,-1],
                    center = TRUE,
                    scale = TRUE)
#saveRDS(data.pca,"/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/pca.features.rds")
#data.pca = readRDS("/data/local/kwreczy/HOT_DATA/generic_ML_for_kmer_hot/pca.features.rds")

#' The plot method returns a plot of the variances (y-axis) associated with the PCs (x-axis).
plot(data.pca, type = "l", main="Variances (y-axis) associated with the PCs (x-axis)")

#' The summary method describe the importance of the PCs.
#summary(data.pca)

data.cols = ifelse(data[,1]>=.99, "red", "black")
data.orgs = rep(names(orgs), len.perc.orgs)
data.shapes = rep(1:4, len.perc.orgs) #1=hg19;2=mm9;2=dm3;3=ce10
stopifnot(length(data.cols) == length(data.shapes))
print( table(data.cols) )
hot.cold.indx = lapply(1:4, function(i){
  hot.indx = intersect(which(data.shapes==i), which(data[,1]>=.99))
  cold.indx = intersect(which(data.shapes==i), which(data[,1]<.5))
  list(hot.indx=hot.indx, cold.indx=cold.indx)
})
names(hot.cold.indx) <- names(orgs)

#' Assign colors to hOT and COL regions and shapes to type of organism
names(orgs)
org.colors.hot = c("yellow", "orange", "blue", "green")
org.colors.cold = c( rgb(153/255,153/255,0), rgb(153/255,76/255,0),
                     rgb(0,0,153/255),rgb(0,102/255,0))
#pie(rep(25,4),  col = org.colors.cold)

colors.pca=rep("black", length(data[,1]))
for(i in 1:4){
  #print(names(orgs)[i])
  colors.pca[ hot.cold.indx[[i]]$hot.indx ]=org.colors.hot[i]
  colors.pca[ hot.cold.indx[[i]]$cold.indx ]=org.colors.cold[i]
}

#' Plot pairs of first 4 principal components
pairs(data.pca$x[, 1:4], col=colors.pca)

#' Plot two first principal components
plot(data.pca$x[, 1], data.pca$x[, 2],
     col=colors.pca,
     #pch=shapes.org,
     main = "HOT and COLD regions", xlab = "PC1", ylab = "PC2")
legend("topright",
       c( paste0("HOT ",names(orgs)), paste0("COLD ",names(orgs))  ),
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       col=c(org.colors.hot, org.colors.cold),
       cex=.8,
       inset = .02)

#' Plot first three principal components
require(scatterplot3d)
scatterplot3d(data.pca$x[, 1:3],
              color = colors.pca, 
              angle = 55,
              main="HOT and COLD regions")
legend("topright",
       c( paste0("HOT ",names(orgs)), paste0("COLD ",names(orgs))  ),
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       col=c(org.colors.hot, org.colors.cold),
       cex=.8)

#' Plot two first principal components using smoothScatter
.jets<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
par(mfrow=c(1,1))
smoothScatter(data.pca$x[,1:2],
              #colramp = .jets,
              main="PCA: projected locations\nfor hg19, mm9, dm3 and ce10\n HOT and COLD regions",
              col = colors.pca)
points(data.pca$x[,1:2],pch=20,
       col=add.alpha(colors.pca, alpha=0.1) )
legend("topright",
       c(paste0(names(orgs), c(" HOT")),paste0(names(orgs), c(" COLD"))),
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       col=c(org.colors.hot, org.colors.cold),
       pt.cex = .5,
       cex=.5)

par(mfrow=c(2,2))
for(i in 1:4){
  
  org=names(orgs)[i]
  
  smoothScatter(data.pca$x[,1:2],
                #colramp = .jets,
                main=org)
  # hot
  points(data.pca$x[hot.cold.indx[[i]][[1]], 1:2],pch=20,
         col=org.colors.hot[i] )
  # cold
  points(data.pca$x[hot.cold.indx[[i]][[2]],1:2],pch=20,
         col=org.colors.cold[i] )
  
  legend("topright",
         c("HOT","COLD"),
         lty=c(1,1), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),
         col=c(org.colors.hot[i], org.colors.cold[i]),
         pt.cex = .5,
         cex=.5)
}

