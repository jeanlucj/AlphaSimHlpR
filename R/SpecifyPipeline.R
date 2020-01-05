#' specifyPipeline function
#'
#' function specify the product pipeline. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the pipeline parameters. bsp is short for breeding sheme parameters
#' @return A list containing objects that specify the product pipeline. This list will determine the number of lists in the records object
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
  
  # Vector of number of checks used in each stage
  # Checks are replicated the same as experimental entries
  nChks <- nEntries / 20

  # Vector of number of reps used in each stage
  nReps <- c(1, 1, 2, 2, 3, 3)
  
  # Vector of error variances estimated from historical data 
  # 200 for SDN is a guess
  errVars <- c(200, 146, 82, 40, 37, 37)
  names(nEntries) <- names(nChks) <- names(nReps) <- names(errVars) <- stageNames
  
  # How many cycles to keep records
  nCyclesToKeepRecords=7
  
  bsp <- c(bsp, mget(ls()))
  return(bsp)
}

#' specifyPopulation function
#'
#' Function to specify the species and population characteristics. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @return A list containing objects that specify the species and population characteristics.
#' 
#' @details Call this function before beginning the simulation
#' 
#' @examples
#' bsp <- specifyPopulation()
#' 
#' @export
specifyPopulation <- function(bsp=NULL){
  # Species characteristics
  # Number of chromosomes
  nChr <- 10
  
  # Population characteristics
  # Number of founders
  nFounders <- 600
  # Number of segregating sites per chromosome
  segSites <- 200
  # Number of QTL per chromosome
  nQTL <- 50
  # Number of observed SNP per chromosome
  nSNP <- 50
  # Initial genetic variance
  genVar <- 40
  # Mean and variance of dominance degree
  meanDD <- 0.3; varDD <- 0.01
  
  bsp <- c(bsp, mget(ls()))
  return(bsp)
}
