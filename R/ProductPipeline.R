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
  records <- with(bsp,{
    curYr <- length(records[[1]])
    for (stage in 1:nStages){
      sourcePop <- last(records[[stage]])
      if (stage==1){ # Stage 1: F1 progeny population: random selection use pop
        entries <- selectInd(sourcePop, nInd=nEntries[stage], use="rand", simParam=SP)
      } else{ # Stage > 1: sort the matrix and make population of best
        idBest <- sourcePop$id[order(sourcePop$pheno, decreasing=T)[1:nEntries[stage]]]
        entries <- records[[1]][curYr + 1 - stage][idBest]
      }
      entries <- setPheno(entries, varE=errVars[stage], reps=nReps[stage]*nLocs[stage], simParam=SP)
      phenoRec <- phenoRecFromPop(entries, bsp, stage)
      records[[stage+1]] <- c(records[[stage+1]], list(phenoRec))
    }
    
    # Remove old records if needed
    if (length(records[[2]]) > nCyclesToKeepRecords) records <- removeOldestCyc(records)
    records
  })#END with bsp
  
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
  records <- with(bsp,{
    curYr <- length(records[[1]])
    phenoDF <- framePhenoRec(records)
    # selPipeAdv has to be given in bsp
    selCrit <- selPipeAdv(phenoDF)
    for (stage in 1:nStages){
      sourcePop <- last(records[[stage]])
      if (stage == 1){ # Stage 1 different: no phenotypes but full Pop-class
        idBest <- sourcePop@id
      } else{
        idBest <- order(selCrit[sourcePop$id,], decreasing=T)[1:nEntries[stage]]
        idBest <- sourcePop$id[idBest]
      }
      entries <- records[[1]][curYr + 1 - stage][idBest]
      entries <- setPheno(entries, varE=errVars[stage], reps=nReps[stage]*nLocs[stage], simParam=SP)
      phenoRec <- phenoRecFromPop(entries, bsp, stage)
      # If provided, add checks to the population
      if(!is.null(checks) & nChks[stage] > 0){
        chkPheno <- setPheno(checks[1:nChks[stage]], varE=errVars[stage], reps=chkReps[stage]*nLocs[stage])
        chkRec <- phenoRecFromPop(chkPheno, bsp, stage, checks=T)
        phenoRec <- bind_rows(phenoRec, chkRec)
      }
      records[[stage+1]] <- c(records[[stage+1]], list(phenoRec))
    }#END 1:nStages

    # Remove old records if needed
    if (length(records[[2]]) > nCyclesToKeepRecords) records <- removeOldestCyc(records)
    records
  })#END with bsp
  
  return(records)
}
