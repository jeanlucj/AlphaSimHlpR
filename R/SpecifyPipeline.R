#' specifyPipeline function
#'
#' function specify the product pipeline. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the pipeline parameters
#' @return A list containing objects that specify the product pipeline. This list will also determine the number of lists in the records object
#' 
#' @details Call this function before beginning the simulation
#' 
#' @examples
#' bsp <- specifyPipeline()
#' 
#' @export
specifyPipeline <- function(bsp=NULL){
  # Number of stages in the product pipeline
  nStages <- 6
  # This is the number of stages simulated and retained in the records
  stageNames <- c("SDN", "CET", "PYT", "AYT", "UY1", "UY2")
  
  # Number of parents in the crossing nursery
  nParents <- 30
  # Number of crosses entering the pipeline
  nCrosses <- 100
  # Number of progeny per cross
  nProgeny <- 30
 
  # Vector of number of number of entries in each stage
  nEntries <- c(nCrosses*nProgeny, 1000, 300, 60, 40, 40)
  names(nEntries) <- trialTypeNames
  
  # Vector of number of checks used in each stage
  # Checks are replicated the same as experimental entries
  nChks <- c(nCrosses*nProgeny, 1000, 300, 60, 40, 40) / 20
  names(nChks) <- trialTypeNames
  max_nChk <- max(nChks)
  
  # Vector of number of reps used in each stage
  nReps <- c(1, 1, 2, 2, 3, 3)
  names(nReps) <- trialTypeNames
  
  # Vector of error variances estimated from historical data 
  # 200 for SDN is a guess
  errVars <- c(200, 146, 82, 40, 37, 37)
  names(errVars) <- trialTypeNames
  
}