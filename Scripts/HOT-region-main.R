

  # /*****************************************************************************
  #  HOT-region-main_script.R
  # 
  #  Published with:
  #  "HOT or not: redefining the origin of high-occupancy target regions".
  #  Authors: Katarzyna Wreczycka, Vedran Franke, Bora Uyar, Ricardo Wurmus, Altuna Akalin
  #
  # ******************************************************************************/  
  

###############################################################################
### Define HOT regions in mouse (mm9), human (hg19), 
### worm (ce10) and fly (dm3)
###############################################################################

./Define_HOT_regions/hot.organism.R # Define HOT regions
./Define_HOT_regions/schema_figure.R # Figure 1A
./Define_HOT_regions/peak_rank_HOT.R # Suppl. Figure 1A
./Define_HOT_regions/TFsonHOT.R # Suppl. Figure 1B

###############################################################################
### Assign HOT regions to genes and expression analysis
###############################################################################

#' Expression profiles of HOT region genes
./Expresssion_on_HOT/checkExpression.R # Suppl. Figure 2A
./Expresssion_on_HOT/GeneExpression.R # Figure 1D

#. GO Biological Processes associated with HOT regions
./GREAT_output/GOterms.R # Suppl. Figure 2B

###############################################################################
### Open chromatin analysis on HOT regions
###############################################################################

./OpenChromatin/DnaseqonHOT.R # Figure 1E

###############################################################################
### Data processing and visualization of murine KO ChIP-seq (mm9) and 
### human DRIP/RDIP-seq (hg19) samples
###############################################################################

#' Plot heatmaps and line plots of KO Chip-seq samples of positively and negatively 
#' associated samples with HOT region score
../Heatmaps/KO_heatmaps.R # Figure 3

#' Analysis of DRIP-seq and RDIP-seq samples
../DRIPseq_G4ChIPseqMEthylation_onHOT/Figure4.R # Figure 4 A,B,C

#' Analysis of G4-quadruplexes
../DRIPseq_G4ChIPseqMEthylation_onHOT/Figure4.R # Figure 4D

###############################################################################
### Sequence analyses of HOT regions and elastic net construction
###############################################################################

library(glmnet)
source("../generic_ML/hot_functions.R")
source("../generic_ML/MLreport_cpgi_functions.R")
../Generic_ML/MLreport_cpgi.R # Figure 2A, Supplement. Figure 3A
../Generic_ML/pca_features.R # Figure 2B


###############################################################################
### Methylation on HOT regions
###############################################################################

../DRIPseq_G4ChIPseqMEthylation_onHOT/Figure4.R # Figure 4 E,F





