#' populationImprovement function
#'
#' function to (do something)
#'
#' @param records [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
populationImprovement <- function(records){
  # Make one population to select out of
  bigPop <- mergePops(records[[2]])
  
  # Select parents among all individuals
  parents <- bigPop[selectParents(records, nParents)]
  records[[1]] <- list(randCross(parents, nCrosses=nCrosses, nProgeny=nProgeny, ignoreGender=T))
  return(records)
}
#' selectParents function
#'
#' function to (do something)
#'
#' @param records [value]
#' @param nToSelect [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
selectParents <- function(records, nToSelect){
  allBLUPs <- iidPhenoEval(framePhenoRec(records))
  nToSelect <- min(nToSelect, nrow(allBLUPs))
  return(rownames(allBLUPs)[order(allBLUPs, decreasing=T)][1:nToSelect])
}
