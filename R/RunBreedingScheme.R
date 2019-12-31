#' runBreedingScheme function
#'
#' function to (do something)
#'
#' @param replication [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
runBreedingScheme <- function(replication){
  
  cat("******", replication, "\n")
  initList <- initializeSimParam()
  founderPop <- initList[[1]]
  SP <<- initList[[2]]
  parents <- newPop(founderPop)
  checks <- selectInd(parents, nInd=max_nChk, use="rand")
  checks <- setPheno(checks)
  
  records <- fillPipeline(parents, checks)
  
  for (cycle in 1:nCycles){
    if (cycle > nCyclesToKeepRecords) records <- removeOldestCyc(records)
    cat(cycle, " ")
    # Run the product pipeline componenent
    records <- productPipeline(records, checks=checks, selectFunc=selectAdvance)
    
    # Run the population improvement componenent
    records <- populationImprovement(records)
  }
  cat("\n")
  recordsLast <<- records
  recordMeans <- sapply(records[-1], function(popList) sapply(popList, meanG))
}
