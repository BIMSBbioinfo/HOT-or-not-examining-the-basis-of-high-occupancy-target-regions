require(GenomicRanges)
require(Biostrings)

#' get OE ratio and GC content for a given set of DNAstrings
GCinfo<-function(bs,wins)
{
  if(class(bs)=="BSgenome"){
    str.set=getSeq(bs,wins)
  }else if(class(bs)=="DNAStringSet"){
    str.set=bs
  }else{
    str.set=getSeq(bs,wins)
  }
  
  di.mat=dinucleotideFrequency(  str.set )
  a.mat =alphabetFrequency( str.set ,baseOnly=TRUE )
  #gcsk=(a.mat[,3]-a.mat[,2])/(a.mat[,3] + a.mat[,2]) #it's GC skew
  exp=(a.mat[,2]*a.mat[,3])/width(str.set )
  obs=di.mat[,"CG"]
  OE=obs/exp
  
  cbind(#gcsk=gcsk,
        #gc=100*(a.mat[,2]+a.mat[,3])/rowSums(a.mat[,1:4])
        #obs=obs/width(str.set),
        oe=OE   )
}

#' GC skew function
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

ATskew<-function(str){
  x=letterFrequencyInSlidingView(str, view.width=1, letters=c("A","T"), OR="|", as.prob=FALSE)
  #sk=(x[,2]-x[,1])/(x[,2]+x[,1])
  #max(sk)-min(sk)
  G=-1*x[,1]+x[,2]
  G=cumsum(G)
  cbind(diff.skew=max(G)-min(G), max.skew=max(abs(G)) )
}

#' function compares oligonucleotide frequencies between two sets of
#' DNA string sets
nucFreqCompare<-function(str.set1,str.set2,width,
                         plot=TRUE,p.cut=0.0001,return.all=FALSE){
  freqs1=oligonucleotideFrequency(str.set1,width)
  freqs2=oligonucleotideFrequency(str.set2,width)
  
  result=data.frame(oligo=character(ncol(freqs1)),
                    p.val=numeric(ncol(freqs1)),logfc=numeric(ncol(freqs1)),
                    stringsAsFactors = FALSE)
  for(i in 1:ncol(freqs1)){
    
    result[i,1]= colnames(freqs1)[i]
    
    result[i,2:3]=c(wilcox.test(freqs1[,i],freqs2[,i])$p.value,
                    log2( (mean(freqs1[,i])+1)/(mean(freqs2[,i])+1))
    )
  }
  
  res2=result[result$p.val<p.cut & result$logfc>0 ,]
  
  if(plot){
    
    ncol=ceiling(sqrt(nrow(res2)))
    nrow=floor(sqrt(nrow(res2)))
    
    layout(matrix(1:(nrow*ncol),ncol=ncol,nrow=nrow))
    par(mar = c(1,3,1,1))
    for(i in which(result$p.val<p.cut & result$logfc>0 )){
      boxplot(freqs1[,i],freqs2[,i],main=result[i,1],outline=FALSE)
    }
  }
  
  if(return.all){
    result 
  }else{
    res2
  }
}

# after scaling in some cell of training data
# NA, NAN etc might show up. I just remove
# columns that contain at least one NA, NAN, etc
is.finite.matrix <- function(obj){
  apply(obj, 2, FUN=function(x) all(is.finite(x)))
}
# after scaling some columns/cells have NA values
# I just remove such columns. Mostly this is Unkown
# genomic repeated element
remove.cols.with.NA <- function(data.list){
  col.remove = c() # which columns will be removed from all dataset from list
  for(i in 1:length(data.list)){
    datalist.i = scale(data.list[[i]]) # scaling
    a = is.finite.matrix( datalist.i)
    which.col.rem = setdiff(1:ncol(datalist.i), which(a))
    col.remove = c(col.remove, which.col.rem)
  }
  col.remove = unique(col.remove)
  if(length(col.remove)==0) return(data.list)

  print("After scaling sampled training datasets some values are not finite (NA, NAN etc).")
  print(paste0("Removing columns: ", colnames(data.list[[1]])[col.remove] ))
  for(i in 1:length(data.list)){
    data.list[[i]] <- data.list[[i]][, -c(col.remove)]
  }
  return(data.list)
}

# create multiple models to remedy class imbalance
# the data.list should contain samples of the 
# class that has more data points
mult_glmnet<-function(data.list,y.list,alpha){
  require(glmnet)
  # initialize output data structure
  models=list()
  importance=list()
  
  data.list = remove.cols.with.NA(data.list)

  for(i in 1:length(data.list)){
    
    # train the model with scaling before the modeling
    # so that coefficients have the same scale
    cv.glmn=cv.glmnet(x=scale(data.list[[i]]),
                      standardize=F,
                      y=y.list[[i]],
                      family="binomial",type.measure="auc",alpha=alpha)
    
    importance[[i]]=coef(cv.glmn)[,1][-1]
    
    # create the final model to be saved
    cv.glmn=cv.glmnet(x= data.list[[i]],
                      y=y.list[[i]],
                      family="binomial",
                      type.measure="auc",alpha=alpha)
    
    models[[i]]=cv.glmn
  }
  importance=do.call("cbind",importance)
  
  # plot the means of coefficients as importance
  plot(rowMeans(importance),pch="")
  text(1:length(rowMeans(importance)),
       rowMeans(importance),rownames(importance) )
  return(list(models=models,imp=importance))
}


# divide data set by
#' @param data data to be sampled
#' @param perc percentiles per row of the data 
#' @param target.perc data between these percentiles are the "yes" class
#' in a classification problem
#' @param downsample.perc data between these percentiles is downsamples 
sampleFromQuantiles<-function(data,perc,
                              target.perc=c(1,0.99),
                              downsample.perc=c(0.8,0.5)){
  target=data[perc>target.perc[2] & perc<= target.perc[1],]
  ctrl  =data[perc>downsample.perc[2] & perc<= downsample.perc[1],]
  if( nrow(target) < nrow(ctrl) ){
    ctrl=ctrl[sample(1:nrow(ctrl),nrow(target)),]
  }
  
  res=rbind(target,ctrl)
  y=c( rep(1,nrow(target)),rep(0,nrow(target))  )
  cbind(y,res)
}




#' The purpouse of this function is to cope with that our
#' classifier unintentionally trains itself to predict CpGi if the training set is not balanced.
#'
#' divide data set by
#' @param data data to be sampled
#' @param perc percentiles per row of the data 
#' @param target.perc data between these percentiles are the "yes" class
#' in a classification problem
#' @param downsample.perc data between these percentiles is downsamples 
sampleFromQuantiles.CpGi.based.sampling<-function(data,perc,
                                                  target.perc=c(1,0.99),
                                                  downsample.perc=c(0.8, 0),
                                                  cpgi, hot){
  require(GenomicRanges)
  require(rtracklayer)
  
  target=data[perc>target.perc[2] & perc<= target.perc[1],]
  ctrl=data[perc>downsample.perc[2] & perc<= downsample.perc[1],]
  
  hot.target = hot[perc>target.perc[2] & perc<= target.perc[1],]  #3276
  hot.ctrl = hot[perc>downsample.perc[2] & perc<= downsample.perc[1],] 
  
  hot.target.cpgi = subsetByOverlaps(hot.target, cpgi) #2525
  perc.overlp.cgi = length(hot.target.cpgi) / length(hot.target) *100
  perc.nooverlp.cgi = 100 - perc.overlp.cgi
  
  
  ## because actually we take the same number of rows of ctrl as for target, 
  ## then only number of windows that should overlap with cpgi is enough
  
  if( nrow(target) < nrow(ctrl) ){
    # adding improvement for cpgi ratio in target and control windows,
    # so then ratio of windows that overlap cpgi islands in control will be the same as in target
    hot.ctrl$id = 1:length(hot.ctrl)
  
    hot.ctrl.cpgi = subsetByOverlaps(hot.ctrl, cpgi)
    hot.ctrl.noncpgi = hot.ctrl[setdiff(hot.ctrl$id, hot.ctrl.cpgi$id),]
    
    
    if( length(hot.ctrl.cpgi) < length( hot.target.cpgi) ){
      #this is the case for mm9, because there is less control regions on cpgi than
      # target regions on cpgi
      # even thought there is a lot control regions, apparently there is not so much cpgi
      hot.target.cpgi.1 = hot.target.cpgi
      hot.target.cpgi = hot.target.cpgi[ sample(1:length(hot.target.cpgi), length(hot.ctrl.cpgi) ), ]
      number_of_ctrl_noncpgi = length(hot.target) - length(hot.target.cpgi)
      
      sample_ctrl_cpgi_ids = sample(hot.ctrl.cpgi$id, length(hot.target.cpgi))
      sample_ctrl_noncpgi_ids = sample(hot.ctrl.noncpgi$id, number_of_ctrl_noncpgi)
      
      ctrl = ctrl[c(sample_ctrl_cpgi_ids, sample_ctrl_noncpgi_ids ),]
      stopifnot(nrow(ctrl) == length(hot.target))
      
    } else{
      # normal case, like before
      
      # size of target and control will be the same
      number_of_ctrl_cpgi = (perc.overlp.cgi * length(hot.target)) / 100
      number_of_ctrl_cpgi=as.numeric(as.character(number_of_ctrl_cpgi)) #mm9 2525
      stopifnot( all.equal(number_of_ctrl_cpgi, length(hot.target.cpgi)) )
      number_of_ctrl_noncpgi = length(hot.target) - number_of_ctrl_cpgi #mm9  751
      number_of_ctrl_noncpgi = as.numeric(as.character(number_of_ctrl_noncpgi))
      
      sample_cpgi_ids = sample(hot.ctrl.cpgi$id, number_of_ctrl_cpgi)
      sample_noncpgi_ids = sample(hot.ctrl.noncpgi$id, number_of_ctrl_noncpgi) #750
      
      ctrl = ctrl[c(sample_cpgi_ids, sample_noncpgi_ids ),]
      stopifnot(nrow(ctrl) == length(hot.target))
      
    }

  }
  res=rbind(target,ctrl)
  y=c( rep(1,nrow(target)),rep(0,nrow(target))  )
  cbind(y,res)
}

#' predict using a set of models
predictFromSetMod<-function(model.list,newx){
  
  row.has.na <- apply(newx, 1, function(x){any(is.na(x))})
  #if(sum(row.has.na) > 0){
  #  print(paste0("Number of rows with some NAs =", sum(row.has.na)))
  #  print("Removing these rows from newx in the predictFromSetMod function")
  #  newx = newx[-which(row.has.na),]  
  #}
  
  require(glmnet)
  pred=matrix(0,ncol=length(model.list)+1,nrow=nrow(newx))
  colnames(pred)=c(paste0("model",1:length(model.list)),"modelAv")
  for(i in 1:length(model.list)){
    pred[,i]=predict(model.list[[i]],data.matrix(newx),
                     type="response")
  }
  pred[,length(model.list)+1]=rowMeans(pred[,1:length(model.list)])
  pred
}


get.BS.genomeSpecies<-function(species.tag){
  if(species.tag=="mm9"){
    BS.genomeSpecies="BSgenome.Mmusculus.UCSC.mm9"
  } else if (species.tag=="hg19"){
    BS.genomeSpecies="BSgenome.Hsapiens.UCSC.hg19"
  } else if (species.tag=="dm3"){
    BS.genomeSpecies="BSgenome.Dmelanogaster.UCSC.dm3"
  } else if (species.tag=="ce10"){
    BS.genomeSpecies="BSgenome.Celegans.UCSC.ce10"
  }
  BS.genomeSpecies
}
