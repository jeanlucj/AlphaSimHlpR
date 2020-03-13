#' removeOldestCyc function
#'
#' function to remove records of the oldest cycles still in \code{records}. Useful to avoid accumulating too much data which slows simulations down and makes them bulky
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @return A \code{records} object with the first population of each list removed.
#' @details \code{records} is a list of lists. This function deletes the first object of each list, excluding the F1 list.
#' 
#' @examples
#' records <- removeOldestCyc(records)
#' 
#' @export
removeOldestCyc <- function(records, nCyclesToKeepRecords){
  # Remove the phenotypic records that are older
  for (i in 2:length(records)){
    nCycStage <- length(records[[i]])
    if (nCycStage > nCyclesToKeepRecords){
      records[[i]] <- records[[i]][-(1:(nCycStage-nCyclesToKeepRecords))]
    }
  }
  # List the id of the individuals remaining
  allID <- NULL
  for (i in 1:length(records[[2]])) allID <- c(allID, records[[2]][[i]]$id)
  for (i in 3:length(records)) allID <- c(allID, records[[i]][[1]]$id)
  allID <- unique(allID)
  allID <- allID[order(as.integer(allID))]
  records[[1]] <- records[[1]][allID]
  return(records)
}
