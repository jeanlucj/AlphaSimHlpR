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
#' @param percentRanges Numeric matrix (nStages + 1 rows and 2 columns) of percentage budget allocation to crossing (F1) and all of the stages.  If there is a stage that is genotyped, the genotyping cost is added to that stage
#' @param startCycle Integer the start cycle from which to measure gain. The end cycle will just be the last cycle
#' @param tolerance Numerical difference between min amd max percentage budgets for all stages
#' @param initializeFunc Function to initialize the breeding program
#' @param baseFile Filename if you want to have progress saved by batch
#' @param productPipeline Function to advance the product pipeline by one generation
#' @param populationImprovement Function to improve the breeding population and select parents to initiate the next cycle of the breeding scheme
#' @param bsp  A list of breeding scheme parameters. It contains pipeline parameters: nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, stageNames, and nCyclesToKeepRecords. It contains species / population parameters: nChr, segSites, nQTL, genVar, meanDD, varDD, nSNP
#' @param randomSeed Integer seed for random number generator
#' @param nCores Integer number of cores to use for parallel simulation
#' @return Numeric matix with all simulations budget allocations, gen mean change, gen std dev change, total cost.
#' 
#' @details A wrapper to repeatedly simulate a scheme with different budget allocations to find optimal allocations
#' 
#' @examples
#' 
#' @export
optimizeByLOESS <- function(batchSize, nByPareto=round(batchSize*0.7), targetBudget, percentRanges, startCycle, tolerance, baseFile=NULL, maxNumBatches=10, initializeFunc, productPipeline, populationImprovement, bsp, randomSeed=1234, nCores=1){
  require(parallel)

  if (length(randomSeed) == batchSize * maxNumBatches){
    randSeeds <- randomSeed
  } else{
    set.seed(randomSeed)
    randSeeds <- round(runif(batchSize * maxNumBatches, min=-1e9, max=1e9))
    saveRDS(randSeeds, file=paste0(baseFile, "randSeeds.rds"))
  }
  
  # Guardrails
  nByPareto <- min(batchSize, nByPareto)
  
  ### Define functions
  # Run the breeding scheme and return the relevant information
  runOneRep <- function(replication, percentRanges, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp, seed=NULL){
    on.exit(expr=saveRDS(list(replication, percentRanges, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp, seed), file="~/runOneRep.rds"))
    if (!is.null(seed)) set.seed(seed)
    bsp$budgetSamplingDone <- FALSE
    while (!bsp$budgetSamplingDone){
      bsp <- sampleEntryNumbers(bsp, targetBudget, percentRanges)
      if (!bsp$budgetSamplingDone){
        # Sampling failed, so shift some budget from later stages to earlier stages
        for (i in nrow(percentRanges):2){
          for (j in 1:2){
            if (percentRanges[i, j] > percentRanges[i-1, j]){
              percentRanges[i, j] <- percentRanges[i, j] - 1
              percentRanges[i-1, j] <- percentRanges[i-1, j] + 1
            }
          }
        }
      }#END make sure proper sampling of budgets was done
    }# Carry on
    rbsOut <- runBreedingScheme(replication=replication, nCycles=bsp$nCyclesToRun, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, bsp=bsp)
    
    on.exit()
    return(list(bsp=bsp, stageOutputs=rbsOut$records$stageOutputs))
  }#END runOneRep
  
  # run a simulation in the vicinity of a previous simulation
  repeatSim <- function(parmRow, replication, radius=0.02, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp, seed=NULL){
    on.exit(expr=saveRDS(list(parmRow, replication, radius, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp, seed), file="~/repeatSim.rds"))
    budg <- parmRow %>% dplyr::select(contains("budget"))
    percentRanges <- t(sapply(unlist(budg), function(prc) c(max(0, prc - radius), min(1, prc + radius))))
    rorOut <- runOneRep(replication, percentRanges, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp, seed)
    return(rorOut)
    on.exit()
  }
  
  # Pull interesting parameters from from the stageOutputs
  getParmsResponse <- function(oneSim, startCycle){
    on.exit(expr=saveRDS(list(oneSim, startCycle), file="~/getParmsResp.rds"))
    bsp <- oneSim$bsp
    parms <- c(budget=bsp$budgetPercentages, nProgeny=bsp$nProgeny, bsp$nEntries, totCost=bsp$totalCosts)
    so <- oneSim$stageOutputs
    resp <- (dplyr::filter(so, stage=="F1" & year==bsp$nCyclesToRun) %>% dplyr::select(genValMean, genValSD)) - (dplyr::filter(so, stage=="F1" & year==startCycle) %>% dplyr::select(genValMean, genValSD))
    on.exit()
    return(unlist(c(parms, resp)))
  }
  
  ### Run optimization
  # Repeatedly simulate the scheme to identify budget allocations optimizing gain
  allBatches <- toRepeat <- tibble()
  allPR <- c(unlist(percentRanges), nSimClose=NA, bestGain=NA, bestSE=NA)
  batchesDone <- 0
  toleranceMet <- FALSE
  while (batchesDone < maxNumBatches & !toleranceMet){
    strtRep <- nrow(allBatches)
    # Repeat a batch of simulations
    if (nrow(toRepeat) > 0){
      saveRDS(toRepeat, file=paste0(baseFile, "toRepeat.rds"))
      #      repeatBatch <- lapply(1:nrow(toRepeat), function(i) repeatSim(toRepeat[i,], strtRep+i, radius=0.04, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, targetBudget=targetBudget, bsp=bsp, seed=randSeeds[strtRep+i]))
      repeatBatch <- mclapply(1:nrow(toRepeat), function(i) repeatSim(toRepeat[i,], strtRep+i, radius=0.04, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, targetBudget=targetBudget, bsp=bsp, seed=randSeeds[strtRep+i]), mc.cores=nCores)
    } else repeatBatch <- list()
    
    strtRep <- strtRep + nrow(toRepeat)
    cat("\n", "@@@@@ nrow(toRepeat)", nrow(toRepeat), "\n")
    # Get a new batch of simulations
#    newBatch <- lapply(1:(batchSize - nrow(toRepeat)), function(i) runOneRep(strtRep+i, percentRanges=percentRanges, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, targetBudget=targetBudget, bsp=bsp, seed=randSeeds[strtRep+i]))
    newBatch <- mclapply(1:(batchSize - nrow(toRepeat)), function(i) runOneRep(strtRep+i, percentRanges=percentRanges, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, targetBudget=targetBudget, bsp=bsp, seed=randSeeds[strtRep+i]), mc.cores=nCores)
    
    # Put together with previous simulatinos
    newBatch <- as_tibble(t(sapply(c(repeatBatch, newBatch), getParmsResponse, startCycle=startCycle)), .name_repair="universal")
    allBatches <- allBatches %>% bind_rows(newBatch)
    
    # Non-Parametric LOESS response
    loFormula <- paste0("genValMean ~ ", paste0(bsp$stageNames, collapse=" + "))
    loFM <- loess(loFormula, data=allBatches)
    loPred <- predict(loFM, se=T)
    bestFit <- which.max(loPred$fit)
    bestSE <- loPred$se.fit[bestFit]
    bestClose <- which(max(loPred$fit) - loPred$fit < 2*bestSE)
    percentRanges <- t(apply(allBatches[bestClose, ] %>% dplyr::select(contains("budget")), 2, range))
    
    # choose which have high response and high std. err. of response
    toRepeat <- tibble()
    fitStdErr <- tibble(batchID=1:nrow(allBatches), fit=loPred$fit, se=loPred$se.fit)
    while (nrow(toRepeat) < nByPareto){
      nonDomSim <- returnNonDom(fitStdErr, dir1Low=F, dir2Low=F, var1name="fit", var2name="se")
      rows <- nonDomSim$batchID
      if (nrow(nonDomSim) > nByPareto - nrow(toRepeat)){
        rows <- sample(rows, nByPareto - nrow(toRepeat))
      }
      toRepeat <- toRepeat %>% bind_rows(allBatches[rows,])
      fitStdErr <- fitStdErr %>% dplyr::filter(!(batchID %in% rows))
    }
    
    allPR <- cbind(allPR, c(unlist(percentRanges), nSimClose=length(bestClose), bestGain=max(loPred$fit), bestSE=bestSE))
    
    # Save batches and results
    if (!is.null(baseFile)){
      saveRDS(allBatches, file=paste0(baseFile, "allBatches.rds"))
      saveRDS(allPR, file=paste0(baseFile, "allPercentRanges.rds"))
    }
    
    batchesDone <- batchesDone + 1
    toleranceMet <- all(percentRanges[,2] - percentRanges[,1] < tolerance)
  }#END keep going until maxNumBatches or tolerance

  return(list(allBatches=allBatches, allPercentRanges=allPR))
}
