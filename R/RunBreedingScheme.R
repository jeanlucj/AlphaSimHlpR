#' runBreedingScheme function
#'
#' function to run a two-part strategy breeding scheme. See Gaynor et al. 2017 for the general idea.
#'
#' @param replication Integer replication of running the breeding scheme
#' @param nCycles Integer number of cycles to run the breeding scheme
#' @param initializeFunc Function to initialize the breeding program.
#' @param productPipeline Function to advance the product pipeline by one generation
#' @param populationImprovement Function to improve the breeding population and select parents to initiate the next cycle of the breeding scheme
#' @param bsp  A list of breeding scheme parameters. It contains pipeline parameters: nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, stageNames, and nCyclesToKeepRecords. It contains species / population parameters: nChr, segSites, nQTL, genVar, meanDD, varDD, nSNP
#' @return A \code{records} object containing the phenotypic records retained of the breeding scheme
#' 
#' @details A wrapper to initiate the breeding program then iterate cycles of product pipeline and population improvement
#' 
#' @examples

#' @export
runBreedingScheme <- function(replication=NULL, nCycles=2, initializeFunc, productPipeline, populationImprovement, bsp){
  cat("******", replication, "\n")
  initList <- initializeFunc(bsp)
  SP <- initList$SP
  bsp <- initList$bsp
  records <- initList$records
  
  for (cycle in 1:nCycles){
    cat(cycle, " ")
    records <- productPipeline(records, bsp, SP)
    records <- populationImprovement(records, bsp, SP)
  }
  cat("\n")
  
  # Finalize the stageOutputs
  records <- lastCycStgOut(records, bsp, SP)
  return(list(records=records, bsp=bsp, SP=SP))
}

#' optimizeByLOESS function
#'
#' function to optimize a two-part strategy breeding scheme:
#' 1. Simulate a batch using given percentage ranges
#' 2. Perform LOESS fit to the gains
#' 3. Find budget with best estimated gain
#' 4. Calculate new percentage ranges: any simulation within 2*StdErr of best
#' Go back to 1.
#'
#' @param batchSize Integer number of simulations between LOESS fits
#' @param targetBudget Numeric value that you want the budget adjusted to
#' @param percentRanges Numeric matrix of percentage budget allocation to crossing (F1) and all of the stages.  If there is a stage that is genotyped, the genotyping cost is added to that stage
#' @param startCycle Integer the start cycle from which to measure gain. The end cycle will just be the last cycle
#' @param tolerance Numerical difference between min amd max percentage budgets for all stages
#' @param initializeFunc Function to initialize the breeding program
#' @param baseFile Filename if you want to have progress saved by batch
#' @param productPipeline Function to advance the product pipeline by one generation
#' @param populationImprovement Function to improve the breeding population and select parents to initiate the next cycle of the breeding scheme
#' @param bsp  A list of breeding scheme parameters. It contains pipeline parameters: nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, stageNames, and nCyclesToKeepRecords. It contains species / population parameters: nChr, segSites, nQTL, genVar, meanDD, varDD, nSNP
#' @param nCores Integer number of cores to use for parallel simulation
#' @return Numeric matix with all simulations budget allocations, gen mean change, gen std dev change, total cost.
#' 
#' @details A wrapper to repeatedly simulate a scheme with different budget allocations to find optimal allocations
#' 
#' @examples
#' 
#' @export
optimizeByLOESS <- function(batchSize, targetBudget, percentRanges, startCycle, tolerance, baseFile=NULL, maxNumBatches=10, initializeFunc, productPipeline, populationImprovement, bsp, nCores=1){
  # Run the breeding scheme and return the relevant information
  runOneRep <- function(replication, percentRanges, initializeFunc, productPipeline, populationImprovement, bsp){
    bsp <- sampleEntryNumbers(bsp, targetBudget, percentRanges)
    rbsOut <- runBreedingScheme(replication=replication, nCycles=bsp$nCyclesToRun, initializeFunc=initializeScheme, productPipeline=productPipeline, populationImprovement=popImprov1Cyc, bsp=bsp)
    return(list(bsp=bsp, stageOutputs=rbsOut$records$stageOutputs))
  }
  
  # Pull interesting parameters from from the stageOutputs
  getParmsResponse <- function(oneSim, startCycle){
    bsp <- oneSim$bsp
    parms <- c(budget=bsp$budgetPercentages, nProgeny=bsp$nProgeny, bsp$nEntries, totCost=bsp$totalCosts)
    so <- oneSim$stageOutputs
    resp <- (filter(so, stage=="F1" & year==bsp$nCyclesToRun) %>% dplyr::select(genValMean, genValSD)) - (filter(so, stage=="F1" & year==startCycle) %>% dplyr::select(genValMean, genValSD))
    return(unlist(c(parms, resp)))
  }
  
  allBatches <- tibble()

  allPR <- c(unlist(percentRanges), nSimClose=NA, bestSE=NA)
  batchesDone <- 0
  toleranceMet <- FALSE
  while (batchesDone < maxNumBatches & !toleranceMet){
    # Get a new batch of simulations
    if (nCores == 1){
      newBatch <- lapply(nrow(allBatches) + 1:batchSize, runOneRep, percentRanges=percentRanges, initializeFunc=initializeScheme, productPipeline=productPipeline, populationImprovement=popImprov1Cyc, bsp=bsp)
    } else{
      newBatch <- mclapply(nrow(allBatches) + 1:batchSize, runOneRep, percentRanges=percentRanges, initializeFunc=initializeScheme, productPipeline=productPipeline, populationImprovement=popImprov1Cyc, bsp=bsp, mc.cores=nCores)
    }
    newBatch <- as_tibble(t(sapply(newBatch, getParmsResponse, startCycle=startCycle)), .name_repair="universal")
    allBatches <- allBatches %>% bind_rows(newBatch)
    
    # Non-Parametric Loess response
    loFormula <- paste0("genValMean ~ ", paste0(bsp$stageNames, collapse=" + "))
    loFM <- loess(loFormula, data=allBatches)
    loPred <- predict(loFM, se=T)
    bestFit <- which.max(loPred$fit)
    bestSE <- loPred$se.fit[bestFit]
    bestClose <- which(max(loPred$fit) - loPred$fit < 2*bestSE)
    percentRanges <- apply(allBatches[bestClose, ] %>% dplyr::select(contains("budget")), 2, range)
    allPR <- cbind(allPR, c(unlist(percentRanges), nSimClose=length(bestClose), bestSE=bestSE))
    
    # Save batches and results
    if (!is.null(baseFile)){
      saveRDS(allBatches, file=paste0(baseFile, "allBatches.rds"))
      saveRDS(allPR, file=paste0(baseFile, "allPercentRanges.rds"))
    }
    
    batchesDone <- batchesDone + 1
    toleranceMet <- all(percentRanges[2,] - percentRanges[1,] < tolerance)
  }#END keep going until stop instructions
  return(list(allBatches=allBatches, allPercentRanges=allPR))
}
