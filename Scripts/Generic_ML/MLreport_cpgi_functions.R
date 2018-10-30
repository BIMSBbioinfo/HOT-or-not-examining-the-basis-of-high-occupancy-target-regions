

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

plot.aveimportance = function(org, mods.path){
  for(species.tag in org){
    # load model 
    mods = readRDS(mods.path[[species.tag]])
    importance = mods$imp
    
    plot(rowMeans(importance),pch="", main=species.tag)
    text(1:length(rowMeans(importance)),
         rowMeans(importance),rownames(importance) ) 
  }    
}


plot.predscores.boxplots = function(comb, mods.path, test.path,
                                    data_path="/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/"){
  
  for(i in 1:nrow(comb)){
    
    x = comb[i,][[1]]
    y = comb[i,][[2]]
    
    print(paste0("test set on ",x,", run model on ", y))
    # train on x organism and predict on y organism 
    # (training data of x organisms, model for y organisms)
    data = readRDS(test.path[[x]])
    mods = readRDS(mods.path[[y]])
    #preds=predictFromSetMod(model.list=mods$models,
    #                        newx=data[,-1])
    #saveRDS(preds, paste0(output_path,"preds.trainon", x,"predicton", y,".rds"))
    preds = readRDS(paste0(output_path,"preds.trainon", x,"predicton", y,".rds"))
    
    # read hot windows
    hot.perc=c(1,0.99)
    downsample.perc=c(0.85,0)
    add.chr=FALSE
    if (x=="ce10") add.chr=TRUE
    if (x=="hg19") hot.perc=c(1,0.995)
    BS.genomeSpecies=get.BS.genomeSpecies(x)
    library(BS.genomeSpecies, character.only = TRUE)
    if(x=="hg19"){
      fname=paste0(data_path, "/HOT", x, "_filterNearBy2kb_Uniform_nopol2a_summits.bed")
    }else if(x=="mm9"){
      fname=paste0(data_path, "/HOT", x, "_nopol2a_summits.bed")
    }else{
      fname=paste0(data_path, "/HOT", x, ".bed")
    }
    
    hot=processHOTsummits(fname=fname,
                          extend=1000,
                          BS.genomeSpecies,
                          add.chr=add.chr)
    # plot prediction scores for windows in organism y
    boxplot(preds[,4][hot$perc<downsample.perc[1]],preds[,4][hot$perc>hot.perc[2]],
            ylab="prediction score",names=c("NOT","HOT"),
            main=paste0(y," model predicting ", x,
                        "\nprediction scores for windows"))
  }
}

# Get true class labels
get.true.labels = function(species.tag, 
                           data_path="/data/akalin/Projects/AAkalin_HOTRegions/Data/HOTregions/"){
  hot.perc=c(1,0.99)
  downsample.perc=c(0.85,0)
  add.chr=FALSE
  if (species.tag=="ce10") add.chr=TRUE
  if (species.tag=="hg19") hot.perc=c(1,0.995)
  BS.genomeSpecies=get.BS.genomeSpecies(species.tag)
  
  if(species.tag=="hg19"){
    fname=paste0(data_path, "/HOT", species.tag, "_filterNearBy2kb_Uniform_nopol2a_summits.bed")
    }else if(species.tag=="mm9"){
      fname=paste0(data_path, "/HOT", species.tag, "_nopol2a_summits.bed")
  }else{
    fname=paste0(data_path, "/HOT", species.tag, ".bed")
  }
  hot=processHOTsummits(fname=fname,
                        extend=1000,
                        BS.genomeSpecies,
                        add.chr=add.chr)
  true_labels = rep(0, length(hot)) # 0 for non-hot windows
  true_labels[which(hot$perc>hot.perc[2])] = 1 # 1 for hot windows
  true_labels
}

plot.hist <- function(ma, xlab, ylab, col){
  require(reshape2)
  require(ggplot2)
  require(RColorBrewer)
  
  melted_cormat <- melt(ma)
  g = ggplot(data = melted_cormat, aes(x=X1, y=X2, fill=value)) + 
    #theme_bw()+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=15))+
    geom_tile(aes(fill = value)) +
    geom_text(aes(fill = melted_cormat$value, label = round(melted_cormat$value, 2)), size=3) +
    scale_fill_gradientn(colours=col, name="AUC")+ 
    xlab(xlab) +
    ylab(ylab)
  #g+theme(axis.text=element_text(size=21),
   #       axis.title=element_text(size=21))
  g
}

# org is a character vector with species tags, e.g. hg19, mm9.
calc.AUC.matrix = function(org, mods.path, test.path){
  require(ROCR)
  ma = matrix(0, ncol=length(org), nrow=length(org))
  colnames(ma) = org; rownames(ma) = org
  for(i in 1:nrow(ma)){
    species.tag.x = rownames(ma)[i] # test data
    for(j in 1:ncol(ma)){
      species.tag.y = colnames(ma)[j] # model
      
      print(paste0("models ", species.tag.y))
      print(paste0("test ", species.tag.x))
      
      
      data = readRDS(test.path[[species.tag.x]])
      mods = readRDS(mods.path[[species.tag.y]])
      
      preds=predictFromSetMod(model.list=mods$models,
                              newx=data[,-1])
      if(species.tag.x!=species.tag.y){
        true.labels = get.true.labels(species.tag.x)
        pred <- prediction( preds[,4], true.labels)
        auc.perf = performance(pred, measure = "auc")
        auc.val = auc.perf@y.values[[1]]
      }else{
        # diagonal of matrix
        model = readRDS(mods.path[[species.tag.x]])
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


#' ## Barplots showing normalized per species importance scores of features to 0-100 scale and average across species
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
    theme(axis.text=element_text(size=10),
             axis.title=element_text(size=15,face="plain")) +
    geom_bar(stat="identity", aes(fill = x), colour="steelblue", position = "dodge", fill="steelblue") + 
    coord_flip() +
    ylab(xlab) + xlab(ylab) + ggtitle(title) 
  #+ scale_y_discrete(breaks=c("0","50","100"),
  #                   labels=c("Dose 0.5", "Dose 1", "Dose 2"))
  #g+theme_bw()
  g
}

plot.betas.barplots = function(org, mods.path, topx=10){
  betas_scaled = list()
  for(species.tag in org){
    
    model=readRDS(mods.path[[species.tag]])
    
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
    
    #plot.barplot(df, species.tag, "Relative importance", "Features")
  }
  
  matr = do.call(cbind,
                 lapply(1:length(betas_scaled), function(i) betas_scaled[[i]]$scaled_beta))
  colnames(matr) = org
  rownames(matr) = as.character(betas_scaled[[1]]$x)
  matr = cbind(matr, means=rowMeans(matr))
  
  df = data.frame(x=rownames(matr),
                  y=matr[,which(colnames(matr) == "means")], 
                  rank=calcRank(matr[,which(colnames(matr) == "means")]))
  df20 = head(df[with(df, order(rank)), ], n=topx)
  df20 <- df[rev(rownames(df20)),]
  plot.barplot(df20, paste0("All models (top ",topx," features)"), "Average relative importance", "Features")

}

