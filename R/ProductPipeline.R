#' prodPipeSimp function
#'
#' Simple function to advance a simulated breeding product pipeline forward by one generation. No selection function, no checks. Selection on phenotype
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has new records created by advancing by a generation
#' 
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
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
prodPipeSimp <- function(records, bsp, SP){
  toAdd <- list()
  for (stage in 1:bsp$nStages){
    if (stage==1){ # Stage 1: F1 progeny population: take them at random
      nGenoRec <- nInd(records[[1]])
      nF1 <- bsp$nCrosses * bsp$nProgeny
      indToAdv <- nGenoRec - nF1 + sort(sample(nF1, bsp$nEntries[1]))
    } else{ # Stage > 1: sort the matrix and make population of best
      sourcePop <- last(records[[stage]])
      indToAdv <- sourcePop$id[order(sourcePop$pheno, decreasing=T)[1:bsp$nEntries[stage]]]
    }
    entries <- records[[1]][indToAdv]
    entries <- setPheno(entries, varE=bsp$errVars[stage], reps=bsp$nReps[stage]*bsp$nLocs[stage], simParam=SP)
    phenoRec <- phenoRecFromPop(entries, bsp, stage)
    toAdd <- c(toAdd, list(phenoRec))
  }
  for (stage in 1 + 1:bsp$nStages){
    records[[stage]] <- c(records[[stage]], toAdd[stage-1])
  }
  
  # Remove old records if needed
  if (length(records[[2]]) > bsp$nCyclesToKeepRecords) records <- removeOldestCyc(records, bsp$nCyclesToKeepRecords)
  
  return(records)
}

#' prodPipeFncChk function
#'
#' function to advance a simulated breeding product pipeline forward by one generation. See Gaynor et al. 2017 for the general idea.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has new records created by advancing by a generation
#' 
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- prodPipeFncChk(records, bsp, SP)
#' records <- popImprov1(records, bsp, SP)
#'
#' @export
prodPipeFncChk <- function(records, bsp, SP){
  phenoDF <- framePhenoRec(records)
  # selPipeAdv has to be given in bsp
  selCrit <- bsp$selPipeAdv(phenoDF)
  toAdd <- list()
  for (stage in 1:bsp$nStages){
    if (stage == 1){ # Stage 1 different: no phenotypes but full Pop-class
      nGenoRec <- nInd(records[[1]])
      nF1 <- bsp$nCrosses * bsp$nProgeny # Sample from the most-recent F1s
      indToAdv <- nGenoRec - nF1 + sort(sample(nF1, bsp$nEntries[1]))
    } else{ # 1:bsp$nEntries[stage-1]] keeps only non-checks
      sourcePop <- last(records[[stage]])
      indToAdv <- order(selCrit[sourcePop$id[1:bsp$nEntries[stage-1]]], decreasing=T)[1:bsp$nEntries[stage]]
      indToAdv <- sourcePop$id[indToAdv]
    }
    entries <- records[[1]][indToAdv]
    entries <- setPheno(entries, varE=bsp$errVars[stage], reps=bsp$nReps[stage]*bsp$nLocs[stage], simParam=SP)
    phenoRec <- phenoRecFromPop(entries, bsp, stage)
    # If provided, add checks to the population
    if(!is.null(bsp$checks) & bsp$nChks[stage] > 0){
      chkPheno <- setPheno(bsp$checks[1:bsp$nChks[stage]], varE=bsp$errVars[stage], reps=bsp$chkReps[stage]*bsp$nLocs[stage], simParam=SP)
      chkRec <- phenoRecFromPop(chkPheno, bsp, stage, checks=T)
      phenoRec <- bind_rows(phenoRec, chkRec)
    }
    toAdd <- c(toAdd, list(phenoRec))
  }#END 1:nStages
  for (stage in 1 + 1:bsp$nStages){
    records[[stage]] <- c(records[[stage]], toAdd[stage-1])
  }

  # Remove old records if needed
  if (length(records[[2]]) > bsp$nCyclesToKeepRecords) records <- removeOldestCyc(records, bsp$nCyclesToKeepRecords)
  
  return(records)
}
