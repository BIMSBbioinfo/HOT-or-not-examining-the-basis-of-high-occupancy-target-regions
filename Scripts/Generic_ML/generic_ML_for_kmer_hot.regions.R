#' # Machine-learning models for HOT regions
#' The models are applied on data sets with k-mer statistics


#  Create training data and run Elastic Net for cross-species comparison of feature conservation in the sequence level.
#' @param species.tag a character that specifies species tag, here one of hg19, mm9, dm3 or ce10 
#' @param data_path a character indicating path to directory where are stored HOT region in bed format
#' @param output_path a character indicating path to directory where training datasets and models should be saved
#' @param cpgiislands_path a character indicating path to directory where are locations of CpG islands in bed format
#' @note
#' Output of this function are two files rds files: training dataset and list of three models and average model
#' 

# data_path="/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/"
# output_path="/data/akalin/kasia/projects/HOT_DATA/generic_ML_for_kmer_hot/"
# cpgiislands_path=paste0("/data/akalin/kasia/projects/151001_UCSC_CpGIslands_",
#                         species.tag,".bed.gz")

createTrainingData.runElasticNet = function(species.tag,
                                             data_path="/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/",
                                             output_path="/data/akalin/kasia/projects/HOT_DATA/generic_ML_for_kmer_hot/",
                                             cpgiislands_path=paste0("/data/akalin/kasia/projects/151001_UCSC_CpGIslands_",
                                                                     species.tag,".bed.gz"),
                                             cpgi_sampling=TRUE # hg19, mm9
                                             ){
  
  
  #' ## Set parameters
  source("/home/kwreczy/projects/HOT-regions/Generic_ML/hot_functions.R")
  
  # which genome should be used
  BS.genomeSpecies=get.BS.genomeSpecies(species.tag)
  
  # percentiles to decide which regions are hot
  hot.perc=c(1,0.99)
  if(species.tag=="hg19") hot.perc=c(1,0.995)
  
  # percentiles to decide which regions are not HOT
  downsample.perc=c(0.85,0)
  
  # if a "chr" string should be added when reading the HOT regions from file
  add.chr=FALSE
  if(species.tag=="ce10") add.chr=TRUE
  
  #' ## Process HOT regions to get features  
  #' read HOT regions from given species
  library(GenomicRanges)
  library(BS.genomeSpecies, character.only = TRUE)
  
  # process regions, get them from file and extend them to 1kb
  if(species.tag=="hg19"){
    hot=processHOTsummits(fname=paste0(data_path, "/HOT",species.tag,"_filterNearBy2kb_Uniform_nopol2a_summits.bed"),
                          extend=1000,
                          BS.genomeSpecies,
                          add.chr=add.chr)  
  }else if(species.tag=="mm9"){
    hot=processHOTsummits(fname=paste0(data_path, "/HOT",species.tag,"_nopol2a_summits.bed"),
                          extend=1000,
                          BS.genomeSpecies,
                          add.chr=add.chr) 
  }else{
    hot=processHOTsummits(fname=paste0(data_path, "/HOT",species.tag,".bed"),
                          extend=1000,
                          BS.genomeSpecies,
                          add.chr=add.chr)
  }
  
  print("create test dataset")
  
  #' get GC info
  BS.genomeSpecies=eval(parse(text=BS.genomeSpecies))
  gcsk=GCinfo(BS.genomeSpecies,hot) #only o/e
  
  #' get sequences from the genome
  seqs=getSeq(BS.genomeSpecies,hot)
  
  #' get GC skew
  #skews=t(sapply(seqs,skew))
  tmp = mclapply(1:length(seqs), skew.paral, str=seqs, mc.cores=50)
  skews = do.call(rbind, tmp)
  
  #' get kmers
  tfrqs2m=oligonucleotideFrequency(seqs,width=2)/(width(seqs)-1)
  tfrqs3m=oligonucleotideFrequency(seqs,width=3)/(width(seqs)-2)
  tfrqs4m=oligonucleotideFrequency(seqs,width=4)/(width(seqs)-3)
  
  #' ## create the data for ML
  #' This is a dataframe where the subset will be used for training 
  data=cbind(gcsk,
             skew=skews[,1],
             tfrqs2m,
             tfrqs3m,
             tfrqs4m)

  
  #' save unsampled training data
  perc=hot$perc
  if(species.tag=="hg19"){
    test.data = paste0(output_path,"HOT.kmer.training",species.tag,"_filterNearBy2kb_Uniform_nopol2a_nogcsk.rds")
    saveRDS(cbind(perc,data),
            test.data)
    data_tmp = cbind(perc,data)
    #data_tmp = readRDS(test.data)
  }else if(species.tag=="mm9"){
    test.data = paste0(output_path,"HOT.kmer.training",species.tag,"_nopol2a_nogcsk.rds")
    saveRDS(cbind(perc,data),
            test.data)
    data_tmp = cbind(perc,data)
    #data_tmp = readRDS(test.data)

  }else{
    test.data = paste0(output_path,"HOT.kmer.training",species.tag,"_nogcsk.rds")
    saveRDS(cbind(perc,data),
            test.data)
    data_tmp = cbind(perc,data)
    #data_tmp = readRDS(test.data)
  }
  
  # Some windows have "N" nucleotides and
  # it appears that glmnet cannot handle NA values
  perc = data_tmp[,1]
  data = data_tmp[,-1]
  row.has.na <- apply(data, 1, function(x){any(is.na(x))})
  if(sum(row.has.na) > 0){
    data = data[-which(row.has.na),]  
    perc = perc[-which(row.has.na)]
    hot = hot[-which(row.has.na),]
  }
  
  #' ## Downsample training data
  #' This balances the training data so that when we train
  #' there are equal number of HOT regions vs NOT-HOT regions
  print("create training dataset")
  
  set.seed(117)
  if( (species.tag %in% c("ce10","dm3") ) | ( species.tag %in% c("hg19","mm9") & cpgi_sampling==FALSE) ){
    samp1=sampleFromQuantiles(data,perc,
                              target.perc=hot.perc,
                              downsample.perc=downsample.perc)
    samp2=sampleFromQuantiles(data,perc,
                              target.perc=hot.perc,
                              downsample.perc=downsample.perc)  
    samp3=sampleFromQuantiles(data,perc,
                              target.perc=hot.perc,
                              downsample.perc=downsample.perc)   
  }
  
  #' Due to the fact that O/E ratio is very predictive for mm9 and hg19
  #' let's take the same ratio of windows that overlap cpgi in target as well in control
  if(species.tag %in% c("hg19","mm9") & cpgi_sampling==TRUE){
    # I use import.bed instead of genomation::readBed, because I dont know how to read bed file without strand.
    cpgi <- import.bed(con=cpgiislands_path)
    set.seed(117)
    samp1=sampleFromQuantiles.CpGi.based.sampling(data,perc,
                                                  target.perc=hot.perc,
                                                  downsample.perc=downsample.perc,
                                                  cpgi, hot)
    samp2=sampleFromQuantiles.CpGi.based.sampling(data,perc,
                                                  target.perc=hot.perc,
                                                  downsample.perc=downsample.perc,
                                                  cpgi, hot)
    samp3=sampleFromQuantiles.CpGi.based.sampling(data,perc,
                                                  target.perc=hot.perc,
                                                  downsample.perc=downsample.perc,
                                                  cpgi, hot)
  }
  
  #' run ElasticNet models for three different training sets
  print("run ElasticNet models for three different training sets")
  data.list=list(samp1[,-1],samp2[,-1],samp3[,-1])
  y.list=list(samp1[,1],samp2[,1],samp3[,1])
  
  
  # run models and store the results in a list
  mods=mult_glmnet(data.list,y.list,alpha=0.05)
 
  
  
  #' save models
  if(species.tag=="hg19" & cpgi_sampling==TRUE){
    models.path = paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag,"0.995_filterNearBy2kb_cpgi_Uniform_nopol2a_nogcsk.rds")
    saveRDS(mods,
            models.path)
  }else if(species.tag=="mm9" & cpgi_sampling==TRUE){
    models.path = paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag,"_cpgi_nopol2a_nogcsk.rds")
    saveRDS(mods,
            models.path)
  }else if(species.tag=="hg19" & cpgi_sampling==FALSE){
      models.path = paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag,"0.995_filterNearBy2kb_Uniform_nopol2a_nogcsk.rds")
      saveRDS(mods,
              models.path)
  }else if(species.tag=="mm9" & cpgi_sampling==FALSE){
      models.path = paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag,"_nopol2a_nogcsk.rds")
      saveRDS(mods,
              models.path)
  }else{
    models.path =paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag,"_nogcsk.rds") 
    saveRDS(mods,
            models.path)
  }
  
  
  return(
    list(modelpath=models.path,
         testpath=test.data)
  )
  
}

