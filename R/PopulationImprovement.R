#' populationImprovement function
#'
#' function to improve a simulated breeding population by one cycle. See Gaynor et al. 2017 for the general idea.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param ppp A list of product pipeline parameters.  It contains nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, trialTypeNames, and nCyclesToKeepRecords.  See \code{initializeFunc} for details
#' @return A records object with a new F1 Pop-class object of progeny coming out of a population improvement scheme
#' 
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' 
#' @examples
#' initList <- initializeFunc()
#' SP <- initList$SP
#' records <- productPipeline(initList$records, ppp=initList$ppp, selectFunc=selectAdvance)
#' records <- populationImprovement(records)
#' recordMeans <- mean_records(records)

#' @export
populationImprovement <- function(records, ppp, SP){
  if (is.null(ppp)) ppp <- list(nParents=30, nCrosses=20, nProgeny=10, nStages=2, errVars=c(4, 2), nReps=c(1, 2), nEntries=nCrosses*nProgeny*c(1, 0.5), nChks=c(0, 0), trialTypeNames=c("PYT", "AYT"), nCyclesToKeepRecords=7, checks=NULL)
  
  records <- with(ppp,{

    # Make one population to select out of
    bigPop <- mergePops(records[[2]])
    
    # Select parents among all individuals
    parents <- bigPop[selectParents(records, ppp)]
    records[[1]] <- list(randCross(parents, nCrosses=nCrosses, nProgeny=nProgeny, ignoreGender=T, simParam=SP))
    records
  })
  
  return(records)
}

#' selectParents function
#'
#' function to (do something)
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param ppp A list of product pipeline parameters.  It contains nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, trialTypeNames, and nCyclesToKeepRecords.  See \code{runBreedingScheme} for details
#' @return Character vector of the ids of the selected individuals
#' @details Accesses all individuals in \code{records} to pick the highest ones
#' @examples none
#' @export
selectParents <- function(records, ppp){
  phenoDF <- framePhenoRec(records)
  allBLUPs <- iidPhenoEval(phenoDF, ppp)
  nToSelect <- min(ppp$nParents, nrow(allBLUPs))
  return(rownames(allBLUPs)[order(allBLUPs, decreasing=T)][1:nToSelect])
}
