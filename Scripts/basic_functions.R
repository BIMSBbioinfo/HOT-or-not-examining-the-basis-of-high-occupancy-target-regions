
####################### Processing file formats

gr2df = function(gr){
  
  data.frame(seqnames=seqnames(gr),
             starts=start(gr)-1,
             ends=end(gr),
             names=c(rep(".", length(gr))),
             scores=gr$scores,
             strands=strand(gr))
}

gr2df.narrowPeak = function(gr){
  data.frame(seqnames=as.character(seqnames(gr)),
             starts=start(gr)-1,
             ends=end(gr),
             names=c(rep(".", length(gr))),
             scores=c(rep("0", length(gr))),
             strands=c(rep("*", length(gr))),
             signalValue=c(rep("0", length(gr))),
             pValue=c(rep("0", length(gr))),
             qValue=c(rep("0", length(gr))),
             peak=gr$peak
  )
}
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

bed2gr <- function(file){
  require(readr)
  df <- read_delim(file, col_types="ciicic", delim="\t",col_names = FALSE,
                   locale=locale(decimal_mark = ",", grouping_mark = "_"))
  colnames(df) <- c("chrom","start", "end", "name", "score", "strand")
  g = makeGRangesFromDataFrame(
    df,
    keep.extra.columns=FALSE,
    starts.in.df.are.0based=FALSE,
    ignore.strand=TRUE)
  g$name = df$name
  g$score =df$score
  g
}

my_readTableFast <- function(filename, header=TRUE,skip=0,sep="\t", nrows=200){
  
  filename <- genomation:::compressedAndUrl2temp(filename)
  if(skip==FALSE){
    skip = 0
  }else if(skip=="auto"){
    skip = genomation:::detectUCSCheader(filename)
  }
  
  if(grepl("^.*(.zip)[[:space:]]*$", filename)){
    tab30rows <- read.zip(filename,
                          sep=sep, 
                          skip=skip,
                          nrows=nrows,
                          header=header,
                          stringsAsFactors=FALSE)
  }else{
    tab30rows <- read.table(file=filename, 
                            sep=sep, 
                            skip=skip,
                            nrows=nrows,
                            header=header,
                            stringsAsFactors=FALSE)
  }
  classes  <- sapply(tab30rows, class)
  cl=""
  for(cla in classes){
    if(cla=="character"){
      cl=paste0(cl, "c")
    }else if(cla=="numeric" | cla=="double"){
      cl=paste0(cl, "d")
    }else if(cla=="integer"){
      cl=paste0(cl, "i")
    }else if(cla=="logical")
      cl=paste0(cl, "l")
  }
  
  df <- read_delim(file=filename, 
                   delim=sep, 
                   skip=skip,
                   col_names=header,
                   col_types=cl,
                   locale=locale(grouping_mark = "_")
  )
  #changing default variables names from read_delim (X[0-9]+) to data.frame (V[0-9]+)
  colnames(df) <- gsub("^X(\\d+)$", "V\\1", colnames(df)) 
  return(as.data.frame(df))
}
# require(R.utils);
# require(genomation)
# reassignInPackage("readTableFast", pkgName="genomation", my_readTableFast);


############################################ Processing HOT regions

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
  #' add percentiles
  mcols(hot)$perc=ecdf(hot$score)(hot$score)
  hot = hot[order(-hot$score)]
  
  hot
}


processHOTsummits_ucsc<-function(fname,extend=1000,bs,add.chr=FALSE){
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
  seqinfo(hot)=keepStandardChromosomes(SeqinfoForUCSCGenome(bs))
  #' remove short ones
  hot=hot[width(GenomicRanges::trim(hot)) == extend,]
  #' add percentiles
  mcols(hot)$perc=ecdf(hot$score)(hot$score)
  hot = hot[order(-hot$score)]
  
  hot
}

add.hot.perc.intervals = function(hot){
  hot$x = "0"
  hot[hot$perc > 0.99 &  hot$perc <=1]$x  = "0.99-1"
  hot[hot$perc > 0.9 &  hot$perc <=0.99]$x = "0.9-0.99"
  hot[hot$perc > 0.8 &  hot$perc <=0.9]$x = "0.8-0.9"
  hot[hot$perc > 0.7 &  hot$perc <=0.8]$x = "0.7-0.8"
  hot[hot$perc > 0.6 &  hot$perc <=0.7]$x = "0.6-0.7"
  hot[hot$perc > 0.5 &  hot$perc <=0.6]$x = "0.5-0.6"
  hot[hot$perc > 0.4 &  hot$perc <=0.5]$x = "0.4-0.5"
  hot[hot$perc <=0.4]$x = "0-0.4"
  a = c("0.99-1", "0.9-0.99", "0.8-0.9", "0.7-0.8", "0.6-0.7", "0.5-0.6", "0.4-0.5", "0-0.4")
  list(hot=hot, interv = rev(a))
}


################################################# Calculating log2(IP/control) by using only reads

normlz <- function(signal, regions){
  param <- ScanBamParam(which=regions)
  cnts.sg=countBam(signal, param=param)$records
  
  message("normalizing... ",date())
  command=paste0("/home/aakalin/.guix-profile/bin/samtools idxstats ",signal," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
  sg.mapped=as.numeric(try(system(command,intern=TRUE)))
  
  cnts.sg=1e6*cnts.sg/sg.mapped
  cnts.sg
}

# function for IP/control enrichment over pre-defined regions
enrich<-function(signal,background,regions){
  
  print(signal)
  print(background)
  
  param <- ScanBamParam(which=regions)
  cnts.sg=countBam(signal, param=param)$records
  cnts.bg=countBam(background, param=param)$records
  
  message("normalizing... ",date())
  command=paste0("/home/aakalin/.guix-profile/bin/samtools idxstats ",signal," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
  sg.mapped=as.numeric(try(system(command,intern=TRUE)))
  command2=paste0("/home/aakalin/.guix-profile/bin/samtools idxstats ",background," | awk 'BEGIN {a=0} {a += $3 } END{print a }'")
  bg.mapped=as.numeric(try(system(command2,intern=TRUE)))
  
  cnts.sg=1e6*cnts.sg/sg.mapped
  cnts.bg=1e6*cnts.bg/bg.mapped
  log2((cnts.sg+1)/(cnts.bg+1))
}



