#' removeOldestCyc function
#'
#' function to (do something)
#'
#' @param records [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
removeOldestCyc <- function(records){
  for (i in 2:length(records)) records[[i]] <- records[[i]][-1]
  return(records)
}
