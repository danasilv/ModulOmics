#' Run ModulOmics 
#'
#' Infer cancer driver modules from genetic alterations and expression unique to the cohort, 
#' and PPI file and regulatory network of H_sapiens. Saves (K-1) files named "ModuleOmicsK.txt" to the output library, 
#' containing the sets and scores of each K.
#'
#' @param mutInput location of the binary genetic alterations file 
#' @param zScoreInput location of the zScored expression file
#' @param PPIInput location of the shortest path score for each pair of proteins, computed based on the PPI network
#' @param TRRUSTInput location of the regulatory network
#' @param pathCode location of the library containg funcsTiMExPPI.R, ilpSolver.R and librariesTiMExPPI.R
#' @param pathOutput location to save output
#' @param quantNo quantile to consider active transcription factor, default value = 3
#' @param topResults number of modules to retreive from each cluster, default value = 5
#' @param noClusters Number of clusters to start to stochastic search in order to locate global maximum, default value = 10
#' @param SETS Number of initial seed modules retreived from the ILP, default value = 200
#' @return A matrix of the detected modules and their single and set scores.
#' @export
run_ModulOmics<-function(mutInput, zScoreInput,  PPIInput, TRRUSTInput, pathCode, pathOutput, quantNo = 3, topResults = 5, noClusters = 10, SETS = 200, LOWER_LIMIT = 2, UPPER_LIMIT = 5){

  source(paste(pathCode,"funcsTiMExPPI.R",sep=""))
  source(paste(pathCode,"ilpSolver.R",sep=""))
  source(paste(pathCode,"librariesTiMExPPI.R",sep=""))
  
  corrMethod<-"spearman"
  pvalThresh<-0.05
  ACTIVE_GENE_CUTOFF_UP = 1
  ACTIVE_GENE_CUTOFF_DOWN = -1
  ACTIVE_PERCENTAGE_PATIENTS = 4
  UNIFORM <- TRUE
  
  
  types <- c('ME','PPI','CoReg','CoExp')
  NUM_MEASURES <- length(types)
  
  mutMatrix<-read.table(mutInput,header=TRUE,sep="\t",stringsAsFactors=F, check.names=FALSE)
  matME<-as.matrix(mutMatrix)
  zScoresMatrix<-read.table(zScoreInput, sep="\t", header = TRUE)
  
  # load PPImatrix
  PPImatrix<-read.table(PPIInput,header=TRUE)
  PPImatrix<-PPImatrix[,-1]
  
  # load REGmatrix
  REGmatrix<-read.table(TRRUSTInput,stringsAsFactors = FALSE)
  if (length(which(duplicated(REGmatrix) | duplicated(REGmatrix[nrow(REGmatrix):1, ])[nrow(REGmatrix):1]==1))>0)
    REGmatrix<-REGmatrix[-which(duplicated(REGmatrix) | duplicated(REGmatrix[nrow(REGmatrix):1, ])[nrow(REGmatrix):1]==1),]
  
  
  # Filter to Active TF: active defines as zscore of over abs(1), for more than a quarter of the patients
  numPatients <- nrow(zScoresMatrix)
  active <-  (colSums((zScoresMatrix >= ACTIVE_GENE_CUTOFF_UP) | (zScoresMatrix <= ACTIVE_GENE_CUTOFF_DOWN)))
  activeTF <- names(active[active > (numPatients/ACTIVE_PERCENTAGE_PATIENTS)])
  RegFilterToActive <- data.frame(activeTF, activeTF, rep(0.5,length(activeTF)))
  colnames(RegFilterToActive) <-  colnames(REGmatrix)
  RegFilterToActive2 <-rbind(REGmatrix, RegFilterToActive)
  RegFilterToActive2 <- RegFilterToActive2[-which(duplicated(RegFilterToActive2)),]
  REGmatrix <- RegFilterToActive2[order(RegFilterToActive2[,1]),]
  
  
  # compute pairwise scores
  pairs<-computePairwiseScoresMod(matME,PPImatrix,REGmatrix,zScoresMatrix,pvalThresh,corrMethod,quantNo)
  
  pairsScore<-(pairs$pairsScores$pairsME+pairs$pairsScores$pairsPPI+pairs$pairsScores$pairsCoReg+
                   pairs$pairsScores$pairsCoExp)/NUM_MEASURES
  diag(pairsScore)<-0
  pairsZScore<-(pairs$pairsZScores$ZpairsME+pairs$pairsZScores$ZpairsPPI+
                    pairs$pairsZScores$ZpairsCoReg+pairs$pairsZScores$ZpairsCoExp)/NUM_MEASURES
  diag(pairsZScore)<-0
  
  # create the input graph for ILP
  gZScore<-graph.adjacency(pairsZScore,mode="undirected",weighted = TRUE)
  
  
  # call ILP solver
  results <- list()
  for (K in LOWER_LIMIT:UPPER_LIMIT) {
    
    initialZILP<-ilpSolver(gZScore, K, SETS)
    initialZSets<-initialZILP$sets
    
    # refine the solution
    refinement<-refineClusters(initialZSets,gorder(gZScore),PPImatrix,zScoresMatrix,matME,REGmatrix, pathCode,corrMethod,pvalThresh, quantNo, topResults, noClusters, NUM_MEASURES)
    
    # Write table
    outSets <- matrix(NA,dim(refinement$results_sets)[1],K)  
    for(i in 1:dim(refinement$results_sets)[1]) {
      indexes <- which(refinement$results_sets[i,] == 1)
      outSets[i,]<-colnames(matME)[indexes]  
    }
    
    outResults <- cbind(outSets, refinement$results_score, refinement$results_MEscore, 
                        refinement$results_PPIscore, refinement$results_CoExpscore,refinement$results_CoRegscore)
    

    colnames(outResults) <- c(paste("gene",seq(1,K),sep=""),"ModuleOmics score","ME score","PPI score","CoExpression score","CoRegulation score")
    write.table(outResults, file = sprintf("%s/ModulOmics%d.txt", pathOutput, K),sep = "\t",quote=FALSE,  row.names = F, col.names = T) 

    results[[K-1]]<-refinement

  }
  return (results)
}
