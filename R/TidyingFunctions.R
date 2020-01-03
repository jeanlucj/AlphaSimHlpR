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
removeOldestCyc <- function(records){
  for (i in 2:length(records)) records[[i]] <- records[[i]][-1]
  return(records)
}
