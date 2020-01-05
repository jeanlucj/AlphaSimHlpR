#' runBreedingScheme function
#'
#' function to run a two-part strategy breeding scheme. See Gaynor et al. 2017 for the general idea.
#'
#' @param replication Integer replication of running the breeding scheme
#' @param nCycles Integer number of cycles to run the breeding scheme
#' @param initializeFunc Function to initialize the breeding program.
#' @param productPipeline Function to advance the product pipeline by one generation
#' @param populationImprovement Function to improve the breeding population and select parents to initiate the next cycle of the breeding scheme
#' @param nFounders Integer of the number of founders of the breeding program
#' @param nChr Integer of the number of chromosomes the simulated species
#' @param segSites Integer of the number of segregating sites per chromosome
#' @param nQTL Integer of the number of causal loci per chromosome
#' @param genVar Real of the genetic variance among founders
#' @param meanDD Real of degree of dominance of QTL
#' @param varDD Real of variance of degree of dominance of QTL
#' @param nSNP Integer of number of SNPs observed per chromosome when a population is genotyped
#' @param nParents Integer of the number of parents among which to make crosses. Default is nParents=30
#' @param nCrosses Integer of how many crosses to make among parents. Default is nCrosses=20
#' @param nProgeny Integer of how many progeny to make per cross. Default is nProgeny=10
#' @param nStages Integer of how many trial stages the product pipeline has. Default is nStages=2
#' @param errVars Real vector of length nStages: the error variances of the trials in each stage. Default is errVars=c(4, 2)
#' @param nReps Integer vector of length nStages: the number of reps used in each trial stage. Default is nReps=c(1, 2)
#' @param nEntries Integer vector of length nStages: the number of entries in trials in each stage. Default is nEntries=nCrosses*nProgeny*c(1, 0.5)
#' @param nChks Integer vector of length nStages: the number of checks to include in trials in each stage. Default is nChks=c(0, 0) (no checks)
#' @param stageNames Character vector of length nStages: the names of the trial types for each stage. Default is stageNames=c("PYT", "AYT")
#' @param nCyclesToKeepRecords Integer of the number of breeding cycles over which to keep records. Simulations get slow and bulky if you keep too many cycles. Records older than nCyclesToKeepRecords will be lost completely. Default is nCyclesToKeepRecords=7 (completely arbitrary)
#' @return A \code{records} object containing the phenotypic records retained of the breeding scheme
#' 
#' @details A wrapper to initiate the breeding program then iterate cycles of product pipeline and population improvement
#' 
#' @examples none
#' 
#' @export
runBreedingScheme <- function(replication=NULL, nCycles=2, initializeFunc, productPipeline, populationImprovement, bsp){
  records <- with(bsp, {
    
    cat("******", replication, "\n")
    initList <- initializeFunc(bsp)
    SP <- initList$SP # SP has to go in the Global Environment
    bsp <- initList$bsp
    records <- initList$records
    
    for (cycle in 1:nCycles){
      cat(cycle, " ")
      records <- productPipeline(records, bsp, selectFunc=selectAdvance, SP)
      records <- populationImprovement(records, bsp, selectFunc=selectParIID, SP)
    }
    cat("\n")
    
    records
  })#END with bsp
  return(records)
}
