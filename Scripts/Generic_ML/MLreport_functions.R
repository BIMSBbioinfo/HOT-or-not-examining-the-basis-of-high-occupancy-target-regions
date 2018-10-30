
# x and y are species.tag, e.g. ce10, mm9, hg19, dm3
train.onX.predict.onY = function(test, model, exclude.cols.train=""){
  species.tag.x = test
  species.tag.y = model
  # load test data
  print("load test data")
  print(paste0(output_path, "HOT.kmer.training",species.tag.x,".rds"))
  data = readRDS(paste0(output_path, "HOT.kmer.training",species.tag.x,".rds"))
  
  # I have to remove some columns from train dataset, because 
  # during running models
  # after scaling sampled training datasets some values are not finite
  # and they are removed from models...
  if(exclude.cols.train!="")
    data = data[,-which(colnames(data) %in% exclude.cols.train)]
  
  percwh = which(colnames(data)=="perc") # [tbd] why in some training dataset there is perc column.., its not bad, but annoying.
  if(length(percwh)>0)
    data = data[-percwh,]
  
  # load model
  print("load model")
  print(paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag.y,".rds"))
  mods = readRDS(paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag.y,".rds"))
  # predict
  preds=predictFromSetMod(model.list=mods$models,
                          newx=data)
  preds
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

# Plot means of coefficients as importance for calc. models
plot.aveimportance = function(org, auxfunc){
  for(species.tag in org){
    # load model 
    mods = readRDS(paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag, auxfunc("",species.tag)$model, ".rds"))
    importance = mods$imp
    
    plot(rowMeans(importance),pch="", main=species.tag)
    text(1:length(rowMeans(importance)),
         rowMeans(importance),rownames(importance) ) 
  }    
}

plot.predscores.boxplots = function(comb, auxfunc){
  for(i in 1:nrow(comb)){
    x = comb[i,][[1]]
    y = comb[i,][[2]]
    print(paste0("test set for ",x,", model for ", y))
    # train on x organism and predict on y organism 
    preds = train.onX.predict.onY(test=paste0(x, auxfunc(x,y)$train), model=paste0(y, auxfunc(x,y)$model))
    saveRDS(preds, paste0(output_path,"preds.trainon", 
                          paste0(x, auxfunc(x,y)$train),"predicton", y, auxfunc(x,y)$model,".rds"))
    #preds = readRDS(paste0(output_path,"preds.trainon", 
    #                       paste0(x, auxfunc(x,y)$train),"predicton", y, auxfunc(x,y)$model,".rds"))
    # read hot windows
    hot.perc=c(1,0.99)
    downsample.perc=c(0.85,0)
    add.chr=FALSE
    if (x=="ce10") add.chr=TRUE
    if (x=="hg19") hot.perc=c(1,0.995)
    BS.genomeSpecies=get.BS.genomeSpecies(x)
    library(BS.genomeSpecies, character.only = TRUE)
    if(x=="hg19"){
      hot=processHOTsummits(fname=paste0(data_path, "/HOThg19_filterNearBy2kb_Uniform.bed"),
                            extend=1000,
                            BS.genomeSpecies,
                            add.chr=add.chr)
    }else{
      hot=processHOTsummits(fname=paste0(data_path, "/HOT", x, ".bed"),
                            extend=1000,
                            BS.genomeSpecies,
                            add.chr=add.chr)
    }
    
    # plot prediction scores for windows in organism y
    boxplot(preds[,4][hot$perc<downsample.perc[1]],preds[,4][hot$perc>hot.perc[2]],
            ylab="prediction score",names=c("NOT","HOT"),
            main=paste0(y," model predicting ", x))
  }
}

get.true.labels = function(species.tag, hotsummitpath){
  hot.perc=c(1,0.99)
  downsample.perc=c(0.85,0)
  add.chr=FALSE
  if (species.tag=="ce10") add.chr=TRUE
  if (species.tag=="hg19") hot.perc=c(1,0.995)
  BS.genomeSpecies=get.BS.genomeSpecies(species.tag)
  
  hot=processHOTsummits(fname=hotsummitpath,
                        extend=1000,
                        BS.genomeSpecies,
                        add.chr=add.chr)
  true_labels = rep(0, length(hot)) # 0 for non-hot windows
  true_labels[which(hot$perc>hot.perc[2])] = 1 # 1 for hot windows
  true_labels
}

get.true.labels.old = function(species.tag, hotpath){
  hot.perc=c(1,0.99)
  downsample.perc=c(0.85,0)
  add.chr=FALSE
  if (species.tag=="ce10") add.chr=TRUE
  if (species.tag=="hg19") hot.perc=c(1,0.995)
  BS.genomeSpecies=get.BS.genomeSpecies(species.tag)
  
  fname = paste0(data_path, "/HOT", species.tag , auxfunc(species.tag,"")$train, ".bed")
  #if (species.tag=="hg19") fname =  paste0(data_path, "/HOT", species.tag, ".bed")
  hot=processHOTsummits(fname=fname,
                        extend=1000,
                        BS.genomeSpecies,
                        add.chr=add.chr)
  true_labels = rep(0, length(hot)) # 0 for non-hot windows
  true_labels[which(hot$perc>hot.perc[2])] = 1 # 1 for hot windows
  true_labels
}

plot.matrix <- function(ma, xlab, ylab, col){
  require(reshape2)
  require(ggplot2)
  require(RColorBrewer)
  
  melted_cormat <- melt(ma, varnames = c("model","test.data")) # rows are models, columns test datasets
  
  g = ggplot(data = melted_cormat, aes(x=test.data, y=model, fill=value)) + 
    geom_tile(aes(fill = value)) +
    geom_text(aes(fill = melted_cormat$value, label = round(melted_cormat$value, 2)), size=7) +
    scale_fill_gradientn(colours=col, name="AUC")+ 
    xlab(xlab) +
    ylab(ylab)
  g+theme(axis.text=element_text(size=21),
          axis.title=element_text(size=21))
}


# Heatmap showing AUC of different combinations of models and test data
# org is a character vector with species tags, e.g. hg19, mm9.
# org is a character vector with species tags, e.g. hg19, mm9.
calc.AUC.matrix = function(org, auxfunc, org.test=org, org.model=org){
  require(ROCR)
  ma = matrix(0, ncol=length(org.test), nrow=length(org.model))
  colnames(ma) = org.test; rownames(ma) = org.model
  
  for(i in 1:length(org.model)){
    st.model = rownames(ma)[i] # species tag for model
    for(j in 1:length(org.test)){
      st.test = colnames(ma)[j] # species tag for test data
      
      testsuffix = paste0(st.test, auxfunc(st.test,st.model)$train)
      modelsuffix = paste0(st.model, auxfunc(st.test,st.model)$model)
      
      # train on x organism and predict on y organism
      preds = train.onX.predict.onY(testsuffix,
                                    modelsuffix) 
      saveRDS(preds,
              paste0(output_path,"preds.trainon",testsuffix,"predicton",modelsuffix,".rds"))
      #preds = readRDS(paste0(output_path,"preds.trainon",testsuffix, "predicton",modelsuffix, ".rds"))
      
      if(st.test!=st.model){
        if(st.test=="hg19"){
          hotpath = "/data/akalin/Projects/AAkalin_HOTRegions/Data/HOThg19_filterNearBy2kb_Uniform.bed"
        }else{
          hotpath = paste0(data_path, "/HOT", species.tag , ".bed")
        }
        true.labels = get.true.labels(st.test, hotpath)
        pred <- prediction( preds[,4], true.labels)
        auc.perf = performance(pred, measure = "auc")
        auc.val = auc.perf@y.values[[1]]
      }else{
        # diagonal of matrix
        model=readRDS(paste0(output_path, "HOT.kmer.sampled.ML.models", modelsuffix ,".rds"))
        # http://stats.stackexchange.com/questions/124288/r-glmnet-cross-validated-auc
        # cv.glmnet fits a whole sequence of models, and will report the auc for all of them. 
        # The max of the cvm sequence is the best model's auc.
        # cvm = mean cross-validated error - a vector of length length(lambda)
        auc.val = max(model$models[[3]]$cvm)
      }
      print(paste0("test data=",st.test, ", model=", st.model, ", auc=", auc.val))
      ma[i,j] = auc.val
    }
  }
  ma
}

calc.AUC.matrix.rep.temp = function(org, auxfunc, org.test=org, org.model=org, exclude.cols.test=NA){
  require(ROCR)
  ma = matrix(0, ncol=length(org.test), nrow=length(org.model))
  colnames(ma) = org.test; rownames(ma) = org.model
  
  for(i in 1:length(org.model)){
    st.model = rownames(ma)[i] # species tag for model
    for(j in 1:length(org.test)){
      st.test = colnames(ma)[j] # species tag for test data
      
      testsuffix = paste0(st.test, auxfunc(st.test,st.model)$train)
      modelsuffix = paste0(st.model, auxfunc(st.test,st.model)$model)
      
      species.tag.x = testsuffix
      species.tag.y = modelsuffix
      ### load test data
      print("load test data")
      print(paste0(output_path, "HOT.kmer.training",species.tag.x,".rds"))
      data = readRDS(paste0(output_path, "HOT.kmer.training",species.tag.x,".rds"))
      
      # I have to remove some columns from train dataset, because 
      # during running models
      # after scaling sampled training datasets some values are not finite
      # and they are removed from models...
      #if(!is.na(exclude.cols.test)){
        ect = exclude.cols.test[[st.model]]
        print(ect)
        data = data[,-which(colnames(data) %in% ect)]
      #}
      percwh = which(colnames(data)=="perc") # [tbd] why in some training dataset there is perc column.., its not bad, but annoying.
      if(length(percwh)>0)
        data = data[,-percwh]
      
      ### load model
      print("load model")
      print(paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag.y,".rds"))
      mods = readRDS(paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag.y,".rds"))
      
      
      ### predict
      # model and test data have to have the same set of features, otherwise I get error:
      # Error in cbind2(1, newx) %*% nbeta : 
      #  Cholmod error 'X and/or Y have wrong dimensions' at file ../MatrixOps/cholmod_sdmult.c, line 90
      # that's why I have to add to my test data features/columns that are in model,
      # I add them as columns with only 0s
      train.modelsuffix = paste0(st.model, auxfunc(st.model,st.model)$train)
      print("data.train.model")
      print(paste0(output_path, "HOT.kmer.training",train.modelsuffix,".rds"))
      data.train.model = readRDS(paste0(output_path, "HOT.kmer.training",train.modelsuffix,".rds"))
      percwh = which(colnames(data.train.model)=="perc") # [tbd] why in some training dataset there is perc column.., its not bad, but annoying.
      if(length(percwh)>0){data.train.model = data.train.model[,-percwh]}

      vars = intersect(colnames(data), colnames(data.train.model))
      data2 = data[,vars]
      sdf = setdiff(colnames(data.train.model), colnames(data))
      mat = matrix(0, ncol=length(sdf), nrow=nrow(data2))
      colnames(mat) = sdf
      data2 = cbind(data2, mat)
      
      preds=predictFromSetMod(model.list=mods$models,
                              newx=data2)
      
      saveRDS(preds,
              paste0(output_path,"preds.trainon",testsuffix,"predicton",modelsuffix,".rds"))
      #preds = readRDS(paste0(output_path,"preds.trainon",testsuffix, "predicton",modelsuffix, ".rds"))
      
      if(st.test!=st.model){
        if(st.test=="hg19"){
          hotpath = "/data/akalin/Projects/AAkalin_HOTRegions/Data/HOThg19_filterNearBy2kb_Uniform.bed"
        }else{
          hotpath = paste0(data_path, "/HOT", st.test , ".bed")
        }
        true.labels = get.true.labels(st.test, hotpath)
        pred <- prediction( preds[,4], true.labels)
        auc.perf = performance(pred, measure = "auc")
        auc.val = auc.perf@y.values[[1]]
      }else{
        # diagonal of matrix
        model=readRDS(paste0(output_path, "HOT.kmer.sampled.ML.models", modelsuffix ,".rds"))
        # http://stats.stackexchange.com/questions/124288/r-glmnet-cross-validated-auc
        # cv.glmnet fits a whole sequence of models, and will report the auc for all of them. 
        # The max of the cvm sequence is the best model's auc.
        # cvm = mean cross-validated error - a vector of length length(lambda)
        auc.val = max(model$models[[3]]$cvm)
      }
      print(paste0("test data=",st.test, ", model=", st.model, ", auc=", auc.val))
      ma[i,j] = auc.val
    }
  }
  ma
}

# Heatmap showing AUC of different combinations of models and test data
# org is a character vector with species tags, e.g. hg19, mm9.
# org is a character vector with species tags, e.g. hg19, mm9.
calc.AUC.matrix.temp = function(org, auxfunc, org.test=org, org.model=org, exclude.cols.test=NA){
  require(ROCR)
  ma = matrix(0, ncol=length(org.test), nrow=length(org.model))
  colnames(ma) = org.test; rownames(ma) = org.model
  
  for(i in 1:length(org.model)){
    st.model = rownames(ma)[i] # species tag for model
    for(j in 1:length(org.test)){
      st.test = colnames(ma)[j] # species tag for test data
      
      testsuffix = paste0(st.test, auxfunc(st.test,st.model)$train)
      modelsuffix = paste0(st.model, auxfunc(st.test,st.model)$model)
      
      # train on x organism and predict on y organism
      #preds = train.onX.predict.onY(testsuffix,
      #                              modelsuffix,
      #                              exclude.cols.test)
      species.tag.x = testsuffix
      species.tag.y = modelsuffix
      # load test data
      print("load test data")
      print(paste0(output_path, "HOT.kmer.training",species.tag.x,".rds"))
      data = readRDS(paste0(output_path, "HOT.kmer.training",species.tag.x,".rds"))
      
      # I have to remove some columns from train dataset, because 
      # during running models
      # after scaling sampled training datasets some values are not finite
      # and they are removed from models...
      if(!is.na(exclude.cols.test)){
        ect = exclude.cols.test[[st.model]]
        print(ect)
        data = data[,-which(colnames(data) %in% ect)]
        
      }
      
      percwh = which(colnames(data)=="perc") # [tbd] why in some training dataset there is perc column.., its not bad, but annoying.
      if(length(percwh)>0)
        data = data[,-percwh]
      
      # load model
      print("load model")
      print(paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag.y,".rds"))
      mods = readRDS(paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag.y,".rds"))
      # predict
      preds=predictFromSetMod(model.list=mods$models,
                              newx=data)
      
      saveRDS(preds,
              paste0(output_path,"preds.trainon",testsuffix,"predicton",modelsuffix,".rds"))
      #preds = readRDS(paste0(output_path,"preds.trainon",testsuffix, "predicton",modelsuffix, ".rds"))
      
      if(st.test!=st.model){
        true.labels = get.true.labels(st.test, auxfunc)
        pred <- prediction( preds[,4], true.labels)
        auc.perf = performance(pred, measure = "auc")
        auc.val = auc.perf@y.values[[1]]
      }else{
        # diagonal of matrix
        model=readRDS(paste0(output_path, "HOT.kmer.sampled.ML.models", modelsuffix ,".rds"))
        # http://stats.stackexchange.com/questions/124288/r-glmnet-cross-validated-auc
        # cv.glmnet fits a whole sequence of models, and will report the auc for all of them. 
        # The max of the cvm sequence is the best model's auc.
        # cvm = mean cross-validated error - a vector of length length(lambda)
        auc.val = max(model$models[[3]]$cvm)
      }
      print(paste0("test data=",st.test, ", model=", st.model, ", auc=", auc.val))
      ma[i,j] = auc.val
    }
  }
  ma
}


# old, it's is ok, but should be transpozed before using it
calc.AUC.matrix.old = function(org, auxfunc){
  require(ROCR)
  ma = matrix(0, ncol=length(org), nrow=length(org))
  colnames(ma) = org; rownames(ma) = org
  for(i in 1:nrow(ma)){
    species.tag.x = rownames(ma)[i] # test data
    for(j in 1:ncol(ma)){
      species.tag.y = colnames(ma)[j] # model

      # train on x organism and predict on y organism
      preds = train.onX.predict.onY(paste0(species.tag.x, auxfunc(species.tag.x,"")$train),
                                    paste0(species.tag.y, auxfunc("",species.tag.y)$model)) 
      saveRDS(preds,
              paste0(output_path,"preds.trainon",species.tag.x, auxfunc(species.tag.x,"")$train,
                     "predicton",species.tag.y,auxfunc("",species.tag.y)$model,".rds"))
      preds = readRDS(
              paste0(output_path,"preds.trainon",species.tag.x, auxfunc(species.tag.x,"")$train, 
                     "predicton",species.tag.y,auxfunc("",species.tag.y)$model, ".rds"))
      
      if(species.tag.x!=species.tag.y){
        true.labels = get.true.labels(species.tag.x, auxfunc)
        pred <- prediction( preds[,4], true.labels)
        auc.perf = performance(pred, measure = "auc")
        auc.val = auc.perf@y.values[[1]]
      }else{
        # diagonal of matrix
        model=readRDS(paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag.x, auxfunc("", species.tag.x)$model ,".rds"))
        # http://stats.stackexchange.com/questions/124288/r-glmnet-cross-validated-auc
        # cv.glmnet fits a whole sequence of models, and will report the auc for all of them. 
        # The max of the cvm sequence is the best model's auc.
        # cvm = mean cross-validated error - a vector of length length(lambda)
        auc.val = max(model$models[[3]]$cvm)
      }
      print(paste0("test data=",species.tag.x, ", model=", species.tag.y, ", auc=", auc.val))
      ma[i,j] = auc.val
    }
  }
  ma
}

# ## Barplots showing normalized per species importance scores of features to 0-100 scale and average across species
do01 <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -
                                                 min(x, na.rm=TRUE))}
calcRank <- function(a){
  sorted = sort(a, decreasing = TRUE)
  ranks = 1:length(a)
  match(names(a),names(sorted))
}
#Emulate ggplot2 default color palette
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
#barplot(1:3, col=ggplotColours(n=3))
plot.barplot = function(df, title, xlab, ylab){
  # make V1 an ordered factor
  df$x <- factor(df$x, levels = df$x)
  color = ggplotColours(n=3)[3]
  g = ggplot(df, aes(x, y)) +
    theme_bw(base_size = 20)+
    geom_bar(stat="identity", aes(fill = x), colour=color,position = "dodge", fill=color) + 
    coord_flip() +
    ylab(xlab) + xlab(ylab) + ggtitle(title) 
  #+ scale_y_discrete(breaks=c("0","50","100"),
  #                   labels=c("Dose 0.5", "Dose 1", "Dose 2"))
  plot(g)
}
#barplot(a, names = names(a), 
#        ylab = "", 
#        xlab = "", horiz = TRUE, las = 1,border = NA, col = 1)

#pdf("/home/kwreczy/HOT/scripts/generic_ML/barplots_hg19.0.995_filterNearBy2kb_cpgi.mm9_cpgi.pdf")
plot.betas.barplots = function(org, auxfunc){
  betas_scaled = list()
  for(species.tag in org){
    
    model=readRDS(paste0(output_path, "HOT.kmer.sampled.ML.models",species.tag, auxfunc("",species.tag)$model, ".rds"))
    
    model.imp = rowMeans(model$imp)
    a = do01(abs(model.imp)) * 100
    df = data.frame(x=names(model.imp), 
                    beta=model.imp, 
                    scaled_beta=a,
                    rank=calcRank(a))
    betas_scaled[[species.tag]] = df
    
    df = head(df[with(df, order(rank)), ], n=20)
    df <- df[rev(rownames(df)),]
    df$y = df$scaled_beta
    
    plot.barplot(df, species.tag, "Relative importance", "Features")
  }
  
  matr = do.call(cbind,
                 lapply(1:length(betas_scaled), function(i) betas_scaled[[i]]$scaled_beta))
  colnames(matr) = org
  rownames(matr) = as.character(betas_scaled[[1]]$x)
  matr = cbind(matr, means=rowMeans(matr))
  
  df = data.frame(x=rownames(matr),
                  y=matr[,which(colnames(matr) == "means")], 
                  rank=calcRank(matr[,which(colnames(matr) == "means")]))
  df20 = head(df[with(df, order(rank)), ], n=20)
  df20 <- df[rev(rownames(df20)),]
  plot.barplot(df20, "All models (top 20 features)", "Average relative importance", "Features")
  df10 = head(df[with(df, order(rank)), ], n=10)
  df10 <- df[rev(rownames(df10)),]
  plot.barplot(df10, "All models (top 10 features)", "Average relative importance", "Features")
}












