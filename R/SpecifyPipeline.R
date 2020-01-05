#' specifyPipeline function
#'
#' function specify the product pipeline. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the pipeline parameters. bsp is short for breeding sheme parameters
#' @param big Logical TRUE for big population FALSE for toy for debugging
#' @return A list containing objects that specify the product pipeline. This list will determine the number of lists in the records object
#' 
#' @details Call this function before beginning the simulation
#' 
#' @examples
#' bsp <- specifyPipeline()
#' 
#' @export
specifyPipeline <- function(bsp=NULL, big=FALSE){
  if (!is.null(bsp)) if (exists("big", bsp)) big <- bsp$big
  # Number of stages in the product pipeline
  nStages <- 6
  # This is the number of stages simulated and retained in the records
  stageNames <- c("SDN", "CET", "PYT", "AYT", "UY1", "UY2")
  
  # Number of parents in the crossing nursery
  nParents <- if_else(big, 30, 10)
  # Number of crosses entering the pipeline
  nCrosses <- if_else(big, 100, 20)
  # Number of progeny per cross
  nProgeny <- if_else(big, 30, 10)

  # Vector of number of number of entries in each stage
  if (big){
    nEntries <- c(nCrosses*nProgeny, 1000, 300, 60, 40, 40)
  } else{
    nEntries <- c(nCrosses*nProgeny, 60, 20, 10, 5, 5)
  }
  
  # Vector of number of checks used in each stage
  # Checks are replicated the same as experimental entries
  nChks <- floor(nEntries / 20)

  # Vector of number of reps used in each stage
  nReps <- c(1, 1, 2, 2, 3, 3)
  
  # Vector of error variances estimated from historical data 
  # 200 for SDN is a guess
  errVars <- c(200, 146, 82, 40, 37, 37)
  names(nEntries) <- names(nChks) <- names(nReps) <- names(errVars) <- stageNames
  
  # How many cycles to keep records
  nCyclesToKeepRecords=7
  
  # Function to advance individuals from one stage to the next
  selPipeAdv <- selectAdvIID
  
  bsp <- c(bsp, mget(setdiff(ls(), c("bsp", "big"))))
  return(bsp)
}

#' specifyPopulation function
#'
#' Function to specify the species and population characteristics. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @param big Logical TRUE for big population FALSE for toy for debugging
#' @return A list containing objects that specify the species and population characteristics.
#' 
#' @details Call this function before beginning the simulation
#' 
#' @examples
#' bsp <- specifyPopulation(bsp)
#' 
#' @export
specifyPopulation <- function(bsp=NULL, big=FALSE){
  if (!is.null(bsp)) if (exists("big", bsp)) big <- bsp$big
  # Species characteristics
  # Number of chromosomes
  nChr <- if_else(big, 10, 1)
  
  # Population characteristics
  # Number of founders
  nFounders <- if_else(big, 600, 50)
  # Number of segregating sites per chromosome
  segSites <- if_else(big, 200, 20)
  # Number of QTL per chromosome
  nQTL <- if_else(big, 50, 5)
  # Number of observed SNP per chromosome
  nSNP <- if_else(big, 50, 5)
  # Initial genetic variance
  genVar <- 40
  # Mean and variance of dominance degree
  meanDD <- 0.3; varDD <- 0.01

  bsp <- c(bsp, mget(setdiff(ls(), c("bsp", "big"))))
  return(bsp)
}
