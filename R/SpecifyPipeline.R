#' specifyPipeline function
#'
#' function specify the product pipeline. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the pipeline parameters. bsp is short for breeding sheme parameters
#' @param ctrlFileName The name of the text file with parameter values controling the simulation. Must include the path to the file. If NULL a toy example simulation will be set up
#' @return A list containing objects that specify the product pipeline. This list will determine the number of lists in the records object
#' 
#' @details Call this function before beginning the simulation
#' 
#' @examples
#' bsp <- specifyPipeline()
#' 
#' @export
specifyPipeline <- function(bsp=NULL, ctrlFileName=NULL){
  if (is.null(ctrlFileName)){ # NULL control file: make toy example
    nStages <- 4 # Number of stages in the product pipeline
    stageNames <- c("SDN", "CET", "PYT", "AYT")
    nParents <- 10 # Number of parents in the crossing nursery
    nCrosses <- 20 # Number of crosses entering the pipeline
    nProgeny <- 10 # Number of progeny per cross
    # Number of number of entries in each stage
    nEntries <- c(nCrosses*nProgeny, 60, 20, 10)
    # Number of checks used in each stage
    # Checks are replicated the same as experimental entries
    nChks <- floor(nEntries / 20)
    nReps <- c(1, 1, 2, 2) # Number of reps used in each stage
    # Error variances estimated from historical data 
    # 200 for SDN is a guess
    errVars <- c(200, 146, 82, 40)
    names(nEntries) <- names(nChks) <- names(nReps) <- names(errVars) <- stageNames
    nCyclesToKeepRecords=4 # How many cycles to keep records
    # Function to advance individuals from one stage to the next
    selPipeAdv <- selectAdvIID
    bsp <- c(bsp, mget(setdiff(ls(), "bsp")))
    #END no control file
  } else{
    ctrlParms <- c("nStages", "stageNames", "nParents", "nCrosses", "nProgeny", "nEntries", "nChks", "nReps" ,"errVars", "nCyclesToKeepRecords", "selPipeAdv")
    ctrlParms <- readControlFile(ctrlFileName, ctrlParms)
    names(ctrlParms$nEntries) <- names(ctrlParms$nChks) <- names(ctrlParms$nReps) <- names(ctrlParms$errVars) <- ctrlParms$stageNames
    ctrlParms$selPipeAdv <- get(ctrlParms$selPipeAdv)
    bsp <- c(bsp, ctrlParms)
  }
  return(bsp)
}

#' specifyPopulation function
#'
#' Function to specify the species and population characteristics. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @param ctrlFileName The name of the text file with parameter values specifying the breeding population. Must include the path to the file. If NULL a toy example simulation will be set up
#' @return A list containing objects that specify the species and population characteristics.
#' 
#' @details Call this function before beginning the simulation
#' 
#' @examples
#' bsp <- specifyPopulation(bsp)
#' 
#' @export
specifyPopulation <- function(bsp=NULL, ctrlFileName=NULL){
  if (is.null(ctrlFileName)){ # NULL control file: make toy example
    # Species characteristics
    nChr <- 1 # Number of chromosomes
    # Population characteristics
    nFounders <- 50 # Number of founders
    segSites <- 20 # Number of segregating sites per chromosome
    nQTL <- 5 # Number of QTL per chromosome
    nSNP <- 5 # Number of observed SNP per chromosome
    genVar <- 40 # Initial genetic variance
    meanDD <- 0.3; varDD <- 0.01 # Mean and variance of dominance degree
    bsp <- c(bsp, mget(setdiff(ls(), "bsp")))
    #END no control file
  } else{
    ctrlParms <- c("nChr", "nFounders", "segSites", "nQTL", "nSNP", "genVar", "meanDD", "varDD")
    ctrlParms <- readControlFile(ctrlFileName, ctrlParms)
    bsp <- c(bsp, ctrlParms)
  }
  return(bsp)
}

#' function to read a text control file
#'
#' The text file should be organized as follows
#' 1. Any text after a comment symbol # will be ignored
#' 2. Control parameter names should be on their own line
#' 3. Parameter values should be on the following line. If multiple parameter values are needed they should be separated by white space but on the same line
#' @param fileName The name of the text file to be read. Must include the path to the file
#' @param ctrlParms A string vector with the names of the control parameters that will be searched in the text file
#' @return A named list of the parameter values read from the control file
#' 
#' @details Call this function before beginning the simulation
#' 
#' @examples
#' params <- readControlFile("./inputDir/ctrlFile.txt", c("nStages", "nParents", "nCrosses"))
#' 
#' @export
readControlFile <- function(fileName, ctrlParms){
  ctrlLines <- readLines(fileName)
  ctrlLines <- sapply(ctrlLines, function(st) strsplit(st, "#", fixed=T)[[1]][1])
  findGetParm <- function(parmName){
    parmRow <- grep(parmName, ctrlLines)+1
    parms <- unlist(strsplit(ctrlLines[parmRow], "[[:space:]]"))
    names(parms) <- NULL
    parmsNum <- suppressWarnings(as.numeric(parms))
    if (!any(is.na(parmsNum))) parms <- parmsNum
    return(parms)
  }
  parms <- lapply(ctrlParms, findGetParm)
  names(parms) <- ctrlParms
  return(parms)
}
