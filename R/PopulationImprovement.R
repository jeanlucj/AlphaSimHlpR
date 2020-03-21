#' popImprov1Cyc function
#'
#' Function to improve a simulated breeding population by one cycle. This version takes phenotyped individuals and crosses them to create new F1
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP The AlphaSimR SimParam object
#' @return A records object with a new F1 Pop-class object of progeny coming out of a population improvement scheme
#' 
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- prodPipeSimp(records, bsp, SP)
#' records <- popImprov1(records, bsp, SP)
#' 
#' @export
popImprov1Cyc <- function(records, bsp, SP){
  # Include current year phenotypes for model training?
  trainRec <- records
  if (!bsp$useCurrentPhenoTrain){
    for (stage in 1+1:bsp$nStages){
      trainRec[[stage]] <- trainRec[[stage]][-length(trainRec[[stage]])]
    }
  }
  # Select parents among all individuals
  candidates <- records[[1]]@id
  crit <- bsp$selCritPopImprov(trainRec, candidates, SP)
  parents <- records[[1]][candidates[order(crit, decreasing=T)[1:bsp$nParents]]]
  progeny <- randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreGender=T, simParam=SP)
  records[[1]] <- c(records[[1]], progeny)
  return(records)
}

#' popImprov2Cyc function
#'
#' Function to improve a simulated breeding population by one cycle. This version does two cycles of predicting F1 individuals and making new F1s
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp List of breeding scheme parameters
#' @param SP The AlphaSimR SimParam object
#' @return A records object with the F1 Pop-class object updated with new progeny coming out of a population improvement scheme
#' 
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- prodPipeSimp(records, bsp, SP)
#' records <- popImprov2Cyc(records, bsp, SP)
#' 
#' @export
popImprov2Cyc <- function(records, bsp, SP){
  # Don't include current year (if specified) for the first cycle
  # but do include it for the second cycle
  useCurrentPhenoTrain <- bsp$useCurrentPhenoTrain
  for (cycle in 1:2){
    trainRec <- records
    if (!useCurrentPhenoTrain){
      for (stage in 1+1:bsp$nStages){
        trainRec[[stage]] <- trainRec[[stage]][-length(trainRec[[stage]])]
      }
    }
    candidates <- records[[1]]@id
    crit <- bsp$selCritPopImprov(trainRec, candidates, SP)
    parents <- records[[1]][candidates[order(crit, decreasing=T)[1:bsp$nParents]]]
    progeny <- randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreGender=T, simParam=SP)
    records[[1]] <- c(records[[1]], progeny)
    useCurrentPhenoTrain <- TRUE
  }
  return(records)
}

#' selCritIID function
#'
#' function to select parents among individuals with phenotypes, assuming individual effects are IID
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents
#' @param SP The AlphaSimR SimParam object (not used, here for uniformity)
#' @return An IID BLUP of the trait of the candidates
#' @details Accesses all individuals in \code{records} to pick the highest among candidates. If candidates do not have records, a random sample is returned
#' 
#' @examples
#' allPop <- mergePops(records[[2]])
#' candidates <- allPop@id
#' parents <- allPop[selectParIID(records, candidates, bsp)]
#' 
#' @export
selCritIID <- function(records, candidates, SP){
  phenoDF <- framePhenoRec(records)
  # Candidates don't have phenotypes so return random vector
  if (!any(candidates %in% phenoDF$id)){ 
    crit <- runif(length(candidates))
    names(crit) <- candidates
  } else{
    crit <- iidPhenoEval(phenoDF)
    crit <- crit[candidates]
  }
  return(crit)
}

#' selCritGRM function
#'
#' function to select parents among individuals with phenotypes, assuming individual effects covary according to a GRM
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents
#' @param bsp A list of product pipeline parameters
#' @param SP The AlphaSimR SimParam object (needed to pull SNPs)
#' @return Character vector of the ids of the selected individuals
#' @details Accesses all individuals in \code{records} to pick the highest ones
#' @examples 
#' candidates <- records[[1]][[1]]@id
#' parents <- records[[1]][[1]][selectParGRM(records, candidates, bsp, SP)]
#' 
#' @export
selCritGRM <- function(records, candidates, SP){
  phenoDF <- framePhenoRec(records)
  if (!any(candidates %in% phenoDF$id)){ 
    crit <- runif(length(candidates))
    names(crit) <- candidates
  } else{
    grm <- makeGRM(records, SP)
    # Remove individuals with phenotypes but who no longer have geno records
    # I am not sure this can happen, but it is a safeguard
    phenoDF <- phenoDF[phenoDF$id %in% rownames(grm),]
    crit <- grmPhenoEval(phenoDF, grm)
    crit <- crit[candidates]
  }
  return(crit)
}
