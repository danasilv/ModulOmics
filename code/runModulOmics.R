rm(list=ls())

# Change parameters:

# format: column names are gene followed by either "-Mut" or "-CNA, rownames are patients. Data as downloaded from TCGA
mutInput<-"/home/pace/bnet/danasilv/hitME/upload/SlimInput/mut.txt"
# format: column names are genes, rownames are patients. Data as downloaded from TCGA
zScoreInput<-"/home/pace/bnet/danasilv/hitME/upload/SlimInput/zScore.txt"
pathOutput<-"/home/pace/bnet/danasilv/hitME/upload/output/"
pathCode<-"/home/pace/bnet/danasilv/hitME/upload/code/"
PPIInput<-"/home/pace/bnet/danasilv/hitME/upload/SlimInput/PPI_SP.txt"
TRRUSTInput<-"/home/pace/bnet/danasilv/hitME/upload/SlimInput/TRRUST.txt"

source(paste(pathCode,"callPipeSP.R",sep=""))

# The function reads cancer cohort and computes ModuleOmics scores
results<-run_ModulOmics(mutInput, zScoreInput, PPIInput, TRRUSTInput, pathCode, pathOutput)
  