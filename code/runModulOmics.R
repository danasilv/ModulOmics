# Change parameters if needed (currently assumes that all files are under the current directory):

# format: column names are gene followed by either "-Mut" or "-CNA, rownames are patients. Data as downloaded from TCGA
mutInput <- "SlimInput/mut.txt"
# format: column names are genes, rownames are patients. Data as downloaded from TCGA
zScoreInput <- "SlimInput/zScore.txt"
pathOutput <- "output/"
pathCode <- "code/"
PPIInput <- "SlimInput/PPI_SP.txt"
TRRUSTInput <- "SlimInput/TRRUST.txt"

source(paste(pathCode, "callPipeSP.R", sep = ""))

# The function reads cancer cohort and computes ModuleOmics scores
results <-
  run_ModulOmics(mutInput,
                 zScoreInput,
                 PPIInput,
                 TRRUSTInput,
                 pathCode,
                 pathOutput)
  
