# function to compute the PPI connectivity of a group of genes
computeConnectivity<-function(PPImatrix,currentGroup,permsGroups)
{
  #   PPImatrix: nxn symmetric matrix, where n is the number of mutations (including TFs) in the binary input matrix, in the same order. entry [i,j] is the probability  that i and j are connected
  #    currentGroup: vector of indices, representing the group for which the connectivity is computed
  #    permsGroups: matrix with rows permutations and columns edges, representing all connected components of the currentGroup.
  if (length(currentGroup)!=length(unique(currentGroup)))
    return(0)
  print(currentGroup)
  noNodes<-length(currentGroup)
  
  # Get list of edges in a complete graph of size currentGroup 
  edgeIndices<-t(combn(noNodes,2))
  
  scoresSubgraphs<-apply(permsGroups,1,function(x){
    currentPerm<-x
    scoreToAdd<-1
    #print(currentPerm)
    # To sum subgraph probability - for each edge chosen to the subgraph multiply by p, for each edge which not chosen multiply by (1-p)
    for (i in 1:length(currentPerm))
    {
      #print(currentPerm[i]==1)
      #print(PPImatrix[currentGroup[edgeIndices[i,]][1],currentGroup[edgeIndices[i,]][2]])
      if (currentPerm[i]==0)
        scoreToAdd<-scoreToAdd*(1-PPImatrix[currentGroup[edgeIndices[i,]][1],currentGroup[edgeIndices[i,]][2]])
      if (currentPerm[i]==1)
        scoreToAdd<-scoreToAdd*PPImatrix[currentGroup[edgeIndices[i,]][1],currentGroup[edgeIndices[i,]][2]]
    }
    #print(scoreToAdd)
    return(scoreToAdd)
  })
  scoreGroup<-sum(scoresSubgraphs)
  print(scoreGroup)
  return(scoreGroup)
}

# function to check if a graph is connected
areNodesConnected <- function(graph, data, extra) 
{
  # input: undirected graph and a vector "group" of size(nodes) marked with 1 for all nodes in the group to be assesed and 0 for all nodes not in the group
  # output: boolean of whether the nodes in "group" are connected.
  
  # define the callback function. The function marks on the vector "group" 
  # the current visited node as 0. Once all the nodes in "group" were 
  # visited the search can stop 
  groupSearch[data['vid'] + 1] <<- 0
  sum(groupSearch) == 0
}

# Compute zScores and scale to (0,1)
normalize<-function(scores)
{
  
  # Z-score
  mean_scores <- mean(scores)
  var_scores <- var(scores)

  if (var_scores == 0) {
    var_scores <- 1
    mean_scores <- 0
  }
  
  zScores <- (scores-mean_scores)/sqrt(var_scores)
  
  # Scale to  (0,1)
  if (min(zScores) != max(zScores)) {
    min_zScores <- min(zScores)
    zScores <- zScores - min_zScores
    max_zScores <- max(zScores)
    if (max_zScores != 0) {zScores <- zScores / max_zScores}
  }
  
  return (zScores)
}


# function to refine the clusters:
# Input:
# firstSets - output of ILP heuristics
# NumGenes - total number of genes to consider
# PPImatrix,REGmatrix,zScoresMatrix,matME,pathStart,corrMethod,pvalThresh - as detailed in function computeExactScore4 (used only for calling computeExactScore4)
# Output: the ith set will be held is results_sets[i,], came from cluster results_cluster[i], has a score of results_score[i], 
# which is composed of results_PPIscor(i), results_MEscore(i), results_CoRegscore(i), results_CoExpscore(i):
#results_sets
#results_score
#results_cluster
#results_PPIscore
#results_MEscore
#results_CoRegscore
#results_CoExpscore
#All setes considered (initial sets, and sets makeing it to a representative) had the followinf scores:
#All_PPI_scores
#All_ME_scores
#All_CoReg_scores
#All_CoExpscore
# Inner conctants:
#MAX_TRIALS_FOR_NEW_SET - number of times sampling in the hoipe to find a set not yet considered for current cluster
#MAX_ITR - number of iterations for sampling random gene to improve cluster's top score 
# TOP_RESULTS - how many best results to retreive per cluster
#UNIFORM_RANDOM_GENE_PROBABILITY - sample gene to improve cluster's top set uniformaly, or biased to within the cluster (for best over all sets assign true, for seperating sets to contain divers group of genes assign false)
refineClusters<-function(firstSets, NumGenes,PPImatrix,zScoresMatrix,matME, REGmatrix, pathStart,corrMethod,pvalThresh, quantNo, TOP_RESULTS, NUM_CLUSTERS, NUM_MEASURES, UNIFORM_RANDOM_GENE_PROBABILITY = TRUE)
{

  SETS <- nrow(firstSets)
  REFRESH_ZSCORES <- 1
  MAX_TRIALS_FOR_NEW_SET <- 10
  MAX_ITR <- 20
  
  AllSets <- firstSets # Keep track of all the sets considered for the zscores: the initial set + any new set being chosen for the results
  AllSetsIndex <- seq.int(nrow(AllSets))
  count_NewSets <- 0

  represent <- matrix(0,NUM_CLUSTERS, NumGenes)
  GeneProbabilityComSum <- matrix(0,NUM_CLUSTERS, NumGenes)
  AllSetsIndexRepresent <- rep (0,NUM_CLUSTERS)
  scoreRepresent <- rep (0,NUM_CLUSTERS)
  PPIscoreRepresent <- rep (0,NUM_CLUSTERS)
  MEscoreRepresent <- rep (0,NUM_CLUSTERS)
  CoRegscoreRepresent <- rep (0,NUM_CLUSTERS)
  CoExpscoreRepresent <- rep (0,NUM_CLUSTERS)
  results_sets <- matrix(, nrow =0, ncol = NumGenes)
  results_AllSetsIndex <- matrix(, nrow =0, ncol = 1)
  results_score <- matrix(, nrow =0, ncol = 1)
  results_PPIscore <- matrix(, nrow =0, ncol = 1)
  results_MEscore <- matrix(, nrow =0, ncol = 1)
  results_CoRegscore <- matrix(, nrow =0, ncol = 1)
  results_CoExpscore <- matrix(, nrow =0, ncol = 1)
  results_cluster <- matrix(, nrow =0, ncol = 1)

    # Cluster the ILP seeds
    clusters <- kmeans(firstSets, NUM_CLUSTERS, algorithm="Lloyd")
    # Score each of the sets: Get the 4 different scores for each set, convert the 4 scores to zcores, and average (same scheme as in the pairwise zscoring)  
  SETS_detailed_scores <- apply(firstSets, 1, function(x){
  
    indGroupNow<-which(x > 0)
  
    #computeExactScoreMod<-function(indicesInME, PPImatrix, REGmatrix, zScoresMatrix, matME, pathStart, corrMethod, pvalThresh, quantNo)
    scoreStructNow<-computeExactScoreMod(indGroupNow,PPImatrix,REGmatrix,zScoresMatrix,matME,pathStart,corrMethod,pvalThresh, quantNo)
    scoreNow<-scoreStructNow
    return(scoreNow)
  })
  
  PPI_scores <- rep (0,SETS)
  ME_scores <- rep (0,SETS)
  CoReg_scores <- rep (0,SETS)
  CoExp_scores <- rep (0,SETS)
  
  for (currSet in 1:SETS){
  
    PPI_scores[currSet] <- SETS_detailed_scores[[currSet]]$scorePPI
    ME_scores[currSet] <- SETS_detailed_scores[[currSet]]$scoreME
    CoReg_scores[currSet] <- SETS_detailed_scores[[currSet]]$scoreCoReg
    CoExp_scores[currSet] <- SETS_detailed_scores[[currSet]]$scoreCoEx 
  }
  
  # Compute the averaged zscore for each set
  PPI_zScores<-normalize(PPI_scores)
  ME_zScores<-normalize(ME_scores)
  CoReg_zScores<-normalize(CoReg_scores)
  CoExp_zScores<-normalize(CoExp_scores)
  
  SETS_PPIscores <- PPI_zScores
  SETS_MEscores <- ME_zScores
  SETS_CoRegscores <- CoReg_zScores
  SETS_CoExpscores <- CoExp_zScores
  SETS_scores <- (PPI_zScores + ME_zScores + CoReg_zScores + CoExp_zScores)/4
  
  
  # For each cluster save the highest scoring group and its score 
  for (i in 1:NUM_CLUSTERS) { 
    scoreRepresent[i] <- max(SETS_scores[clusters$cluster == i])
    rep_index <- which(clusters$cluster == i)[which.max(SETS_scores[clusters$cluster == i])]
    represent[i, ] <- firstSets[rep_index,]
    AllSetsIndexRepresent[i] <- rep_index
    PPIscoreRepresent[i] <- SETS_PPIscores[rep_index]
    MEscoreRepresent[i] <- SETS_MEscores[rep_index]
    CoRegscoreRepresent[i] <- SETS_CoRegscores[rep_index]
    CoExpscoreRepresent[i] <- SETS_CoExpscores[rep_index]
  }
  
  # Prepare gene probabilty according to number of appearances in cluster
  if (!UNIFORM_RANDOM_GENE_PROBABILITY) {
    for (i in 1:NUM_CLUSTERS) { 
      if (is.null(dim(firstSets[which(clusters$cluster == i), ]))) {  # There's only ine set in the cluster - no need to colSum
        GeneProbability <- firstSets[which(clusters$cluster == i),]
      } else
      {
        GeneProbability <- colSums(firstSets[which(clusters$cluster == i),])
      }
      GeneProbability <- GeneProbability + 0.1
      SumGeneProbability = sum(GeneProbability)
      GeneProbabilityNormalized <- GeneProbability / SumGeneProbability
      GeneProbabilityComSum[i, ] <- cumsum(GeneProbabilityNormalized)
    }
  }
  # Get TOP RESULTS by iteratively refining each cluster MAX_ITR iterations.
  for (currRound in 1:TOP_RESULTS) {
    for (currC in 1:NUM_CLUSTERS) {
      if (sum(clusters$cluster == currC) > 0) {
        # Check if current cluster merged with previous one - if representative of current cluster is already in the results
        if ((duplicated(rbind(results_sets, represent[currC,]))[nrow(results_sets)+1])) { # Check if was already chosen

        }
        # Check if current cluster has a represntetive (maybe we exhausted representatives)
        else if (any(represent[currC, ] == -1)){
        }
        else
        {
       for (itr in 1:MAX_ITR) { # ToDo: add that if there was no change we also stop?
            representative_set <- which(represent[currC,] > 0)

            # Choose a random gene in the representative group
            replaced_gene = sample (representative_set,1)

            # Choose a random gene in the network
            if (UNIFORM_RANDOM_GENE_PROBABILITY) {
              candidate_gene = sample (1:NumGenes,1)
            } else {
              candidate_gene <- min(which(GeneProbabilityComSum[currC,] >= runif(1, 0, 1)))
            }
            
            # Make sure the candidate gene is not already in the set
            while (represent[currC,][candidate_gene] > 0) {
              candidate_gene = sample (1:NumGenes,1)
            }

            # Replace the gene in the candidate set
            candidate_set <- represent[currC,]
            candidate_set[replaced_gene] = 0
            candidate_set[candidate_gene] = 1
            
            # Make sure the new group is not one that is already assesed (old representatives)
            already_in_results <- duplicated(rbind(results_sets, candidate_set))[nrow(results_sets)+1]

            # If the new group is already in results, Draw new random representative until MAX_TRIALS_FOR_NEW_SET
            trial = 1
            while (already_in_results && trial < MAX_TRIALS_FOR_NEW_SET) {

              # Choose a random gene in the network
              if (UNIFORM_RANDOM_GENE_PROBABILITY) {
                candidate_gene = sample (1:NumGenes,1)
              }
              else {
                candidate_gene <- min(which(GeneProbabilityComSum[currC,] >= runif(1, 0, 1)))
              }
              
              # Make sure the candidate gene is not already in the set
              while (represent[currC,][candidate_gene] > 0) {
                candidate_gene = sample (1:NumGenes,1)
              }
              
              # Replace the gene in the candidate set
              candidate_set <- represent[currC,]
              candidate_set[replaced_gene] = 0
              candidate_set[candidate_gene] = 1
              
              # Make sure the new group is not one that is already assesed (old representatives)
              already_in_results <- duplicated(rbind(results_sets, candidate_set))[nrow(results_sets)+1]
              
              trial = trial + 1
            }
            
            # Got a new set as representative candidate
            if (!already_in_results) {

              # Score the new set, if the score benefits the representative set, make the exchange
              # The new represntative of the cluster did not necessarily come from the ILP.
              new_detailed_score = computeExactScoreMod(which(round(candidate_set,4)==1),PPImatrix,REGmatrix,zScoresMatrix,matME,pathStart,corrMethod,pvalThresh, quantNo)
              
              PPI_zScore <- normalizeOneScore(PPI_scores, new_detailed_score$scorePPI)
              ME_zScore <- normalizeOneScore(ME_scores, new_detailed_score$scoreME)
              CoReg_zScore <- normalizeOneScore(CoReg_scores, new_detailed_score$scoreCoReg)
              CoExp_zScore <- normalizeOneScore(CoExp_scores, new_detailed_score$scoreCoExp)
              
              new_score <- (PPI_zScore+ ME_zScore + CoReg_zScore + CoExp_zScore)/4
              
              # New set is better than the representative - exchane with the representative
              if (new_score > scoreRepresent[currC]) {

                # Save set, score, and index in AllSets
                scoreRepresent[currC] <- new_score
                represent[currC,] <- candidate_set
                row.is.a.match <- apply(AllSets, 1, identical, candidate_set)
                
                # If set exists save its index
                if (any(row.is.a.match)) {
                  AllSetsIndexRepresent[currC] <- which(row.is.a.match)
                }
                else {
                  AllSetsIndexRepresent[currC] = -1 
                }
              }
            }
          }

          # The best set of the cluster is saved, a new representative will be chosen for next iteration
          results_sets <- rbind(results_sets, represent[currC,])
          results_score <- rbind(results_score, scoreRepresent[currC])
          results_PPIscore <- rbind(results_PPIscore, PPIscoreRepresent[currC])
          results_MEscore <- rbind(results_MEscore, MEscoreRepresent[currC])
          results_CoRegscore <- rbind(results_CoRegscore, CoRegscoreRepresent[currC])
          results_CoExpscore <- rbind(results_CoExpscore, CoExpscoreRepresent[currC])
          results_cluster <- rbind(results_cluster, currC)
          
          # If the representative is an existing set 
          if (duplicated(rbind(AllSets,  represent[currC,]))[nrow(AllSets)+1]){
            results_AllSetsIndex <- rbind(results_AllSetsIndex, AllSetsIndexRepresent[currC]) 
          }
          #if the best set is a new set (could happen only if the sample gene created a set that was not in the first sets) - 
          # then the zscores must be recalculated:
          # we save the new scores, and every REFRESH_ZSCORES we recalculate the mean and var and all the scores of all the sets
          else { # Check if was already chosen

            # Append the new set and the new detailed scores
            PPI_scores <- append(PPI_scores,new_detailed_score$scorePPI)
            ME_scores <- append(ME_scores,new_detailed_score$scoreME)
            CoReg_scores <- append(CoReg_scores,new_detailed_score$scoreCoReg)
            CoExp_scores <- append(CoExp_scores,new_detailed_score$scoreCoEx)
            
            AllSets <- rbind(AllSets, represent[currC,])
            results_AllSetsIndex <- rbind(results_AllSetsIndex, dim(AllSets)[1]) 
            
            count_NewSets <- count_NewSets + 1

            # every REFRESH_ZSCORES we recalculate the mean and var and all the scores of all the sets
            if (count_NewSets == REFRESH_ZSCORES) {

              count_NewSets <- 0 
              
              # Update zScores
              PPI_zScores <- normalize(PPI_scores)
              ME_zScores <- normalize(ME_scores)
              CoReg_zScores <- normalize(CoReg_scores)
              CoExp_zScores <- normalize(CoExp_scores)
              
              # Need to update the results sets, the representative sets, and the initial sets still waiting to be chosen from
              tmp_PPI_zScores <- normalizeOneScore(PPI_scores,results_PPIscore)
              tmp_ME_zScores <- normalizeOneScore(ME_scores, results_MEscore)
              tmp_CoReg_zScores <- normalizeOneScore(CoReg_scores, results_CoRegscore)
              tmp_CoExp_zScores <- normalizeOneScore(CoExp_scores, results_CoExpscore)
              
              results_score <- (tmp_PPI_zScores + tmp_ME_zScores + tmp_CoReg_zScores + tmp_CoExp_zScores)/4
              
              tmp_PPI_zScores <- normalizeOneScore(PPI_scores,SETS_PPIscores)
              tmp_ME_zScores <- normalizeOneScore(ME_scores,SETS_MEscores)
              tmp_CoReg_zScores <- normalizeOneScore(CoReg_scores,SETS_CoRegscores)
              tmp_CoExp_zScores <- normalizeOneScore(CoExp_scores,SETS_CoExpscores)
              
              
              SETS_scores <- (tmp_PPI_zScores + tmp_ME_zScores + tmp_CoReg_zScores + tmp_CoExp_zScores)/4
              
              
              # For each cluster update the reprersentative score -  alternatively we can choose a new representative 
              for (i in 1:NUM_CLUSTERS) { 

                tmp_PPI_zScores <- normalizeOneScore(PPI_scores,PPIscoreRepresent[i])
                tmp_ME_zScores <- normalizeOneScore(ME_scores,MEscoreRepresent[i])
                tmp_CoReg_zScores <- normalizeOneScore(CoReg_scores,CoRegscoreRepresent[i])
                tmp_CoExp_zScores <- normalizeOneScore(CoExp_scores,CoExpscoreRepresent[i])
                
                scoreRepresent[i] <- (tmp_PPI_zScores+tmp_ME_zScores+tmp_CoReg_zScores+tmp_CoExp_zScores) / 4
              }
            }
          }
          
          # Choose the next represenatative.
          
          
          # If the old representative is not chosen, it could be it, if it was then the next in the cluster, unless it was also chosen while walking from another cluster
          found_new_representative = FALSE
          while(!found_new_representative && (sum(clusters$cluster == currC) > 0)) {

            score_new_representative <- max(SETS_scores[clusters$cluster == currC])  # Get current maximum scored sets
            set_index <- which(clusters$cluster == currC)[which.max(SETS_scores[clusters$cluster == currC])]
            # If there are still more sets to choose from
            if (is.null(dim(firstSets[which(clusters$cluster == currC), ]))) {

              found_new_representative = TRUE
              represent[currC, ] <- rep(-1,)
            } 
            else  
            {

              set_new_representative <- if (is.matrix(firstSets)) firstSets[set_index,] else firstSets
              if (!(duplicated(rbind(results_sets, set_new_representative))[nrow(results_sets)+1])) { # Check if was already chosen
                
                # Found the new representative
                PPIscoreRepresent[currC] <- SETS_PPIscores[set_index]
                MEscoreRepresent[currC] <- SETS_MEscores[set_index]
                CoRegscoreRepresent[currC] <- SETS_CoRegscores[set_index]
                CoExpscoreRepresent[currC] <- SETS_CoExpscores[set_index]
                scoreRepresent[currC] <- score_new_representative
                represent[currC, ] <- set_new_representative
                row.is.a.match <- apply(AllSets, 1, identical, set_new_representative)
                
                # If set exists save its index
                if (any(row.is.a.match)) {
                  AllSetsIndexRepresent[currC] <- which(row.is.a.match)
                }
                else {
                  
                  # EROOR: the new representative should come from the sets
                  AllSetsIndexRepresent[currC] = -1 
                }
                
                found_new_representative = TRUE
              }
              # The set is already in the results
              else {

                # Remove the set 
                firstSets <- if (is.matrix(firstSets)) firstSets[-set_index,] else -1
                
                SETS_PPIscores <- SETS_PPIscores[-set_index]
                SETS_MEscores <- SETS_MEscores[-set_index]
                SETS_CoRegscores <- SETS_CoRegscores[-set_index]
                SETS_CoExpscores <- SETS_CoExpscores[-set_index]
                
                SETS_scores <- SETS_scores[-set_index]
                clusters$cluster  <- clusters$cluster[-set_index] 
              }
            }
          }
        }
      }
    }
  }
  resultsStruct<-list()
  resultsStruct$results_sets<-results_sets
  resultsStruct$results_score<-results_score
  resultsStruct$results_cluster<-results_cluster
  resultsStruct$results_PPIscore<-results_PPIscore
  resultsStruct$results_MEscore<-results_MEscore
  resultsStruct$results_CoRegscore<-results_CoRegscore
  resultsStruct$results_CoExpscore<-results_CoExpscore
  resultsStruct$All_PPIscore<- PPI_scores
  resultsStruct$All_MEscore<-ME_scores
  resultsStruct$All_CoRegscore<-CoReg_scores
  resultsStruct$All_CoExpscore<-CoExp_scores
  return(resultsStruct)
}



# function to only use the mutations in a TiMEx dataset
onlyMutations<-function(mat)
{
  keep<-which(sapply(colnames(mat),function(x){strsplit(x,"-")}[[1]][2])=="Mut")
  if (length(keep)>0)
    mat<-mat[,keep]
  return(mat)
}



# function to keep both mutations and CNAs, as long as they don't belong to the same gene
keepBoth<-function(mat)
{
  idxMuts<-which(sapply(colnames(mat),function(x){strsplit(x,"-")}[[1]][2])=="Mut")
  nameMuts<-sapply(colnames(mat),function(x){strsplit(x,"-")}[[1]][1])[idxMuts]
  
  idxCNAs<-which(sapply(colnames(mat),function(x){strsplit(x,"-")}[[1]][2])=="CNA")
  nameCNAs<-sapply(colnames(mat),function(x){strsplit(x,"-")}[[1]][1])[idxCNAs]
  
  commonNames<-intersect(nameMuts,nameCNAs)
  if (length(commonNames)>0)
  {
    toRemove<-idxCNAs[match(commonNames,nameCNAs)]
    mat<-mat[,-toRemove]
  }
}

# modified function (CoExp score) to compute the exact score of a group, as a sum of the 4 scores
# quantNo is the quantile to be kept on the absolute z scores, above which correlation will be computed between genes
computeExactScoreMod<-function(indicesInME, PPImatrix, REGmatrix, zScoresMatrix, matME, pathStart, corrMethod, pvalThresh, quantNo)
{
  
  # find ME score
  cat ("4.1")
  groupToTest<-colnames(matME)[indicesInME]
  MEstruct<-testCliqueAsGroup(indicesInME,matME)
  if (MEstruct$pvalueLRT>=pvalThresh){
    scoreME<-0
  } else {scoreME<-MEstruct$opMu$par[length(groupToTest)+1]}
  cat ("4.2\n")
  # find PPI score
  groupPPI<-unique(unlist(lapply(strsplit(groupToTest,"-"),function(x){return(x[1])})))
  cat ("4.21\n")
  groupPPI<-unique(groupPPI)
  cat ("4.22\n")
  notFoundPPI<-groupPPI[which(is.na(match(groupPPI,colnames(PPImatrix)))==1)]
  cat ("4.23\n")
  if (length(notFoundPPI)>0)
    scorePPI<-0
  else
  {
    cat ("4.23\n")
    noGenesPPI<-length(groupPPI)
    if (noGenesPPI==1)
      scorePPI<-0
    else
    {
      permsGroups<-read.table(paste(pathStart,"/data/connected_subgraphs/",noGenesPPI,".txt",sep=""))
      scorePPI<-computeConnectivity(PPImatrix,match(groupPPI,colnames(PPImatrix)),permsGroups)  
    }
  }

  # find CoReg score
  uniqueTFs<-unique(REGmatrix[,1]) # the unique transcription factors
  uniqueTargets<-unique(REGmatrix[,2]) # the unique targets
  targetsMutated<-intersect(groupPPI,uniqueTargets) # the genes in the group which are also targets
  # for each TF, which of the mutated targets it regulates
  matchTFs<-matrix(0,nrow=length(targetsMutated),ncol=dim(REGmatrix)[1])
  colnames(matchTFs)<-REGmatrix[,1]
  rownames(matchTFs)<-targetsMutated
  # for each mutated target, whether each other transcription is a regulator
  matchesUniqueTFs<-matrix(0,nrow=length(targetsMutated),ncol=length(uniqueTFs))
  colnames(matchesUniqueTFs)<-uniqueTFs
  rownames(matchesUniqueTFs)<-targetsMutated

  if (length(targetsMutated)<=1)
  {
    
    scoreCoReg<-0
    structCoReg<-list()                      
  } else {

    for (i in 1:length(targetsMutated))
      matchTFs[i,which(REGmatrix[,2]==targetsMutated[i])]<-1
    for (i in 1:length(uniqueTFs))
    {

      eqnow<-which(colnames(matchTFs)==uniqueTFs[i])
      if (length(eqnow)>1)
      {
        if (dim(matchTFs)[1]==1)
          matchesUniqueTFs[,i]<-sum(matchTFs[,eqnow])
        else
          matchesUniqueTFs[,i]<-apply(matchTFs[,eqnow],1,sum)
      }
      else
        matchesUniqueTFs[,i]<-matchTFs[,eqnow]
    }
  
    # this matrix should have no element greater than 1
    if (length(which(matchesUniqueTFs>1))>0)
      stop("duplicate entries in the regulatory network")
  
    coregulations<-apply(matchesUniqueTFs,2,sum)
    maxCorregulations<-max(coregulations)
    allMax<-which(coregulations==maxCorregulations)
  
    idxMax<-allMax[1]
    if (maxCorregulations>1) {

      scoreCoReg <- maxCorregulations/ length(indicesInME) 
    }
    else
      scoreCoReg<-0
    structCoReg<-list("coregulations"=coregulations,"maxCorregulations"=maxCorregulations,"allMax"=allMax)
  }
  
  # find CoExp score
  notFoundExp<-groupPPI[which(is.na(match(groupPPI,colnames(zScoresMatrix)))==1)]
  noRealGenes<-length(groupPPI)-length(notFoundExp)
  groupExp<-groupPPI
  if (length(notFoundExp)>0)
  {
    groupExp<-setdiff(groupPPI,notFoundExp)
  }
  if (length(groupExp)>1)
  {
  
    absZScores<-apply(zScoresMatrix,2,function(x){mean(abs(x))})
    quantilesNow<-quantile(absZScores)
    minZScore<-quantilesNow[quantNo]
    
    indicesInCorr<-match(groupExp,colnames(zScoresMatrix))
    pairwiseCorr<-matrix(NA,nrow=length(groupExp),ncol=length(groupExp))
    for (g1 in 1:(length(indicesInCorr)-1))
      for (g2 in (g1+1):length(indicesInCorr))
      {
        if ((absZScores[indicesInCorr[g1]]>=minZScore) && (absZScores[indicesInCorr[g2]]>=minZScore))
          pairwiseCorr[g1,g2]<- abs(cor(zScoresMatrix[,indicesInCorr[g1]],zScoresMatrix[,indicesInCorr[g2]],method=corrMethod))
        else
          pairwiseCorr[g1,g2]<-0
      }

    scoreCoExp<-mean(c(as.vector(pairwiseCorr)[!is.na(as.vector(pairwiseCorr))],rep(0,length(notFoundExp))))
  } else
  {
    scoreCoExp<-0
    pairwiseCorr<-NA
  }
  
  return(list("scoreGroup"=(scorePPI+scoreME+scoreCoReg+scoreCoExp)/4,
              "scorePPI"=scorePPI, "scoreME"=scoreME, "scoreCoReg"=scoreCoReg,
              "scoreCoExp"=scoreCoExp, "structME"=MEstruct, "structCoReg"=structCoReg,
              "pairwiseCorr"=pairwiseCorr))
}

# modified function (CoExp score) to compute all 4 pairwise scores
# quantNo is the quantile to be kept on the absolute z scores, above which correlation will be computed between genes
computePairwiseScoresMod<-function(matME,PPImatrix,REGmatrix,zScoresMatrix,pvalThresh,corrMethod,quantNo)
{
  noGenesME<-dim(matME)[2]
  
  # find ME pairwise score
  structsPairsME<-analyzePairs(matME)
  pairsME<-structsPairsME$muEstSym
  pairsPvalME<-structsPairsME$pvalueLRTCorrectSym$uncorrected
  if (length(which(pairsPvalME>pvalThresh))>0)
    pairsME[which(pairsPvalME>pvalThresh)]<-0
  
  # compute ME z-score
  zPairsME<-matrix(NA,nrow = noGenesME, ncol = noGenesME)
  allME<-pairsME[upper.tri(pairsME,diag=FALSE)]
  zPairsME[upper.tri(zPairsME)]<-normalize(allME)
  #zPairsME[upper.tri(zPairsME)]<-(pairsME[upper.tri(zPairsME)]-mean(allME))/sqrt(var(allME))
  #zPairsME[upper.tri(zPairsME)]<- zPairsME[upper.tri(zPairsME)] - min(zPairsME[upper.tri(zPairsME)])
  #if(max(zPairsME[upper.tri(zPairsME)]) != 0) {zPairsME[upper.tri(zPairsME)]<- zPairsME[upper.tri(zPairsME)] / max(zPairsME[upper.tri(zPairsME)])}
  ind<-lower.tri(zPairsME) 
  zPairsME[ind] <- t(zPairsME)[ind] 
  
  
  # find PPI pariwise score
  pairsPPI<-matrix(nrow=noGenesME,ncol=noGenesME)
  colnames(pairsPPI)<-colnames(pairsME)
  rownames(pairsPPI)<-rownames(pairsME)
  notFoundPPI<-c()
  
  for (gene1 in 1:(noGenesME-1))
  {
    print(paste("current gene PPI:", gene1))
    indexG1PPI<-match(strsplit(colnames(matME)[gene1],"-")[[1]][1],colnames(PPImatrix))
    if (is.na(indexG1PPI))
    {
      notFoundPPI<-c(notFoundPPI,strsplit(colnames(matME)[gene1],"-")[[1]][1])
      pairsPPI[gene1,c((gene1+1):noGenesME)]<-rep(0,(noGenesME-gene1))
    }
    
    else
      for (gene2 in (gene1+1):noGenesME)
      {
        indexG2PPI<-match(strsplit(colnames(matME)[gene2],"-")[[1]][1],colnames(PPImatrix))
        if (is.na(indexG2PPI))
        {
          pairsPPI[gene1,gene2]<-0
          notFoundPPI<-c(notFoundPPI,strsplit(colnames(matME)[gene2],"-")[[1]][1])
        }
        else
          pairsPPI[gene1,gene2]<-PPImatrix[indexG1PPI,indexG2PPI]
      }
  }
  notFoundPPI<-unique(notFoundPPI)
  
  # compute PPI z-score
  zPairsPPI<-matrix(NA,nrow = noGenesME, ncol = noGenesME)
  allPPI<-pairsPPI[upper.tri(pairsPPI,diag=FALSE)]
  zPairsPPI[upper.tri(zPairsPPI)]<-normalize(allPPI)
  #zPairsPPI[upper.tri(zPairsPPI)]<-(pairsPPI[upper.tri(zPairsPPI)]-mean(allPPI))/sqrt(var(allPPI)) # Ask Simona about this line
  #zPairsPPI[upper.tri(zPairsPPI)]<- zPairsPPI[upper.tri(zPairsPPI)] - min(zPairsPPI[upper.tri(zPairsPPI)])
  #if (max(zPairsPPI[upper.tri(zPairsPPI)]) != 0) {zPairsPPI[upper.tri(zPairsPPI)]<- zPairsPPI[upper.tri(zPairsPPI)] / max(zPairsPPI[upper.tri(zPairsPPI)])}
  ind<-lower.tri(zPairsPPI) 
  zPairsPPI[ind] <- t(zPairsPPI)[ind] 
  
  
  # find CoReg pairwise score
  pairsCoReg<-matrix(nrow=noGenesME,ncol=noGenesME)
  colnames(pairsCoReg)<-colnames(pairsME)
  rownames(pairsCoReg)<-rownames(pairsME)
  notFoundCoReg<-c()
  uniqueTargets<-unique(REGmatrix[,2])
  
  for (gene1 in 1:(noGenesME-1))
  {
    print(paste("current gene Regulation:", gene1))
    indexG1CoReg<-match(strsplit(colnames(matME)[gene1],"-")[[1]][1],uniqueTargets)
    if (is.na(indexG1CoReg))
    {
      pairsCoReg[gene1,c((gene1+1):noGenesME)]<-rep(0,(noGenesME-gene1))
      notFoundCoReg<-c(notFoundCoReg,strsplit(colnames(matME)[gene1],"-")[[1]][1])
    }
    else
    {
      gene1CoReg<-uniqueTargets[indexG1CoReg]
      for (gene2 in (gene1+1):noGenesME)
      {
        indexG2CoReg<-match(strsplit(colnames(matME)[gene2],"-")[[1]][1],uniqueTargets)
        if (is.na(indexG2CoReg))
        {
          notFoundCoReg<-c(notFoundCoReg,strsplit(colnames(matME)[gene2],"-")[[1]][1])
          pairsCoReg[gene1,gene2]<-0
        }
        else
        {
          gene2CoReg<-uniqueTargets[indexG2CoReg]
          TFs1<-REGmatrix[which(REGmatrix[,2]==gene1CoReg),1]
          TFs2<-REGmatrix[which(REGmatrix[,2]==gene2CoReg),1]
          commonTFs<-intersect(TFs1,TFs2)
          if (length(commonTFs)>0)
            pairsCoReg[gene1,gene2]<-1
          else 
            pairsCoReg[gene1,gene2]<-0
        }
      }
    }
  }
  notFoundCoReg<-unique(notFoundCoReg)
  
  # compute CoReg z-score
  zPairsCoReg<-matrix(NA,nrow = noGenesME, ncol = noGenesME)
  allCoReg<-pairsCoReg[upper.tri(pairsCoReg,diag=FALSE)]
  zPairsCoReg[upper.tri(zPairsCoReg)]<-normalize(allCoReg)
  #zPairsCoReg[upper.tri(zPairsCoReg)]<-(pairsCoReg[upper.tri(zPairsCoReg)]-mean(allCoReg))/sqrt(var(allCoReg))
  #zPairsCoReg[upper.tri(zPairsCoReg)]<- zPairsCoReg[upper.tri(zPairsCoReg)] - min(zPairsCoReg[upper.tri(zPairsCoReg)])
  #if (max(zPairsCoReg[upper.tri(zPairsCoReg)]) != 0) {zPairsCoReg[upper.tri(zPairsCoReg)]<- zPairsCoReg[upper.tri(zPairsCoReg)] / max(zPairsCoReg[upper.tri(zPairsCoReg)])}
  ind<-lower.tri(zPairsCoReg) 
  zPairsCoReg[ind] <- t(zPairsCoReg)[ind] 
  
  
  # find CoExp pariwise score
  pairsCoExp<-matrix(nrow=noGenesME,ncol=noGenesME)
  colnames(pairsCoExp)<-colnames(pairsME)
  rownames(pairsCoExp)<-rownames(pairsME)
  notFoundCoExp<-c()
  
  absZScores<-apply(zScoresMatrix,2,function(x){mean(abs(x))})
  quantilesNow<-quantile(absZScores)
  minZScore<-quantilesNow[quantNo]
  
  for (gene1 in 1:(noGenesME-1))
  {
    print(paste("current gene Expression:", gene1))
    indexG1CoExp<-match(strsplit(colnames(matME)[gene1],"-")[[1]][1],colnames(zScoresMatrix))
    if (is.na(indexG1CoExp))
    {
      pairsCoExp[gene1,c((gene1+1):noGenesME)]<-rep(0,(noGenesME-gene1))
      notFoundCoExp<-c(notFoundCoExp,strsplit(colnames(matME)[gene1],"-")[[1]][1])
    }
    
    else
      for (gene2 in (gene1+1):noGenesME)
      {
        indexG2CoExp<-match(strsplit(colnames(matME)[gene2],"-")[[1]][1],colnames(zScoresMatrix))
        if (is.na(indexG2CoExp))
        {
          pairsCoExp[gene1,gene2]<-0
          notFoundCoExp<-c(notFoundCoExp,strsplit(colnames(matME)[gene2],"-")[[1]][1])
        }
        else
        {
          if ((absZScores[indexG1CoExp]>=minZScore) && (absZScores[indexG2CoExp]>=minZScore))
            pairsCoExp[gene1,gene2]<- abs(cor(zScoresMatrix[,indexG1CoExp],zScoresMatrix[,indexG2CoExp],method=corrMethod))
          else
            pairsCoExp[gene1,gene2]<-0
        }
      }
  }
  notFoundCoExp<-unique(notFoundCoExp)
  
  # compute CoExp z-score
  zPairsCoExp<-matrix(NA,nrow = noGenesME, ncol = noGenesME)
  allCoExp<-pairsCoExp[upper.tri(pairsCoExp,diag=FALSE)]
  zPairsCoExp[upper.tri(zPairsCoExp)]<-normalize(allCoExp)
  #zPairsCoExp[upper.tri(zPairsCoExp)]<-(pairsCoExp[upper.tri(zPairsCoExp)]-mean(allCoExp))/sqrt(var(allCoExp))
  #zPairsCoExp[upper.tri(zPairsCoExp)]<- zPairsCoExp[upper.tri(zPairsCoExp)] - min(zPairsCoExp[upper.tri(zPairsCoExp)])
  #if (max(zPairsCoExp[upper.tri(zPairsCoExp)]) != 0) {zPairsCoExp[upper.tri(zPairsCoExp)]<- zPairsCoExp[upper.tri(zPairsCoExp)] / max(zPairsCoExp[upper.tri(zPairsCoExp)])}
  ind<-lower.tri(zPairsCoExp) 
  zPairsCoExp[ind] <- t(zPairsCoExp)[ind] 
  
  
  # create return structures
  pairsScores<-list()
  pairsScores$pairsME<-pairsME
  pairsScores$pairsPPI<-pairsPPI
  pairsScores$pairsCoReg<-pairsCoReg
  pairsScores$pairsCoExp<-pairsCoExp
  
  pairsZScores<-list()
  pairsZScores$ZpairsME<-zPairsME
  pairsZScores$ZpairsPPI<-zPairsPPI
  pairsZScores$ZpairsCoReg<-zPairsCoReg
  pairsZScores$ZpairsCoExp<-zPairsCoExp
  
  notFound<-list()
  notFound$notFoundPPI<-notFoundPPI
  notFound$notFoundCoReg<-notFoundCoReg
  notFound$notFoundCoExp<-notFoundCoExp
  
  return(list("pairsScores"=pairsScores,"pairsZScores"=pairsZScores,
              "structsPairsME"=structsPairsME,"notFound"=notFound))   
}

# Compute zScores and scale to (0,1)
normalizeOneScore<-function(scores, new_score)
{
  
  # Z-score
  mean_scores <- mean(scores)
  var_scores <- var(scores)
  if (var_scores == 0) {
    var_scores <- 1
    mean_scores <- 0
  }
  
  zScores <- (scores-mean_scores)/sqrt(var_scores)
  new_score <- (new_score - mean_scores)/sqrt(var_scores)
  
  # Scale to  (0,1)
  if (min(zScores) != max(zScores)) {
    min_zScores <- min(zScores)
    zScores <- zScores - min_zScores
    new_score <- new_score - min_zScores
    max_zScores <- max(zScores)
    if (max_zScores != 0) {
      zScores <- zScores / max_zScores
      new_score <- new_score/max_zScores
    }
  }
  
  return (new_score)
}
