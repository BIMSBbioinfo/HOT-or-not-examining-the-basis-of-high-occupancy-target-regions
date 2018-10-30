### HOT regions caller

#' # get peak summit coverage and 2001bp sliding window density
#' @author Altuna Akalin 
getPeakSummitcov = function(tfs, p=501){
  # get all peaks in one GRanges object
  peaks=unlist(tfs)
  end(peaks)=start(peaks) + peaks$peak
  start(peaks)=end(peaks)
  
  # get peak summit coverage and k bp sliding window density
  pcov=coverage(peaks)
  #pcov$chrM <- NULL
  pcov500=runsum(pcov, k = p,endrule="constant")
  pcov500
}

#' # call peaks in a given rle
#' @author Altuna Akalin 
getPeakRle<-function(my.rle,filterNearBy=0){
  require(GenomicRanges)
  
  # get the peaks of TFBS density
  res=GRanges()
  for(chr in names(my.rle)){
    cat(chr,"\n")
    x = diff(sign(diff(c(-Inf, runValue(my.rle[[chr]]), -Inf)))) # tutaj ucina te wszsytkie powtarzajace sie liczby... 

    # get peak locations
    starts=start(my.rle[[chr]])[x == -2]
    ends=end(my.rle[[chr]])[x == -2]
    mid=round((starts+ends)/2) # get the point nearest to all peaks
    values=runValue(my.rle[[chr]])[x == -2]
    
    sset=GRanges(seqnames=chr,ranges=IRanges(mid,mid),scores=values)
    if(filterNearBy>0){
      
      p.rle= coverage(sset,weight=sset$scores)
      # find the max pos in each window, on and rle of peaks only
      maxs=viewWhichMaxs(Views(p.rle[[chr]],start(sset)-filterNearBy,end(sset)+filterNearBy))
      sset=sset[start(sset) == maxs,]
    }
    
    res=c(res,sset )
  }
  res
}
