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
#' @param nCores Integer number of cores to use for parallel simulation
#' @return Numeric matix with all simulations budget allocations, gen mean change, gen std dev change, total cost.
#' 
#' @details A wrapper to repeatedly simulate a scheme with different budget allocations to find optimal allocations
#' 
#' @examples
#' 
#' @export
optimizeByLOESS <- function(batchSize, nByPareto=round(batchSize*0.7), targetBudget, percentRanges, startCycle, tolerance, baseFile=NULL, maxNumBatches=10, initializeFunc, productPipeline, populationImprovement, bsp, nCores=1){
  require(foreach)
  require(doParallel)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  # Guardrails
  nByPareto <- min(batchSize, nByPareto)
  
  ### Define functions
  # Run the breeding scheme and return the relevant information
  runOneRep <- function(replication, percentRanges, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp){
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
    rbsOut <- runBreedingScheme(replication=replication, nCycles=bsp$nCyclesToRun, initializeFunc=initializeScheme, productPipeline=productPipeline, populationImprovement=popImprov1Cyc, bsp=bsp)
    return(list(bsp=bsp, stageOutputs=rbsOut$records$stageOutputs))
  }#END runOneRep
  
  # run a simulation in the vicinity of a previous simulation
  repeatSim <- function(parmRow, replication, radius=0.02, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp){
    budg <- parmRow %>% dplyr::select(contains("budget"))
    percentRanges <- t(sapply(unlist(budg), function(prc) c(max(0, prc - radius), min(1, prc + radius))))
    runOneRep(replication, percentRanges, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp)
  }
  
  # Pull interesting parameters from from the stageOutputs
  getParmsResponse <- function(oneSim, startCycle){
    bsp <- oneSim$bsp
    parms <- c(budget=bsp$budgetPercentages, nProgeny=bsp$nProgeny, bsp$nEntries, totCost=bsp$totalCosts)
    so <- oneSim$stageOutputs
    resp <- (filter(so, stage=="F1" & year==bsp$nCyclesToRun) %>% dplyr::select(genValMean, genValSD)) - (filter(so, stage=="F1" & year==startCycle) %>% dplyr::select(genValMean, genValSD))
    return(unlist(c(parms, resp)))
  }
  
  ### Run optimization
  # Repeatedly simulate the scheme to identify budget allocations optimizing gain
  allBatches <- tibble()
  allPR <- c(unlist(percentRanges), nSimClose=NA, bestGain=NA, bestSE=NA)
  batchesDone <- 0
  toleranceMet <- FALSE
  while (batchesDone < maxNumBatches & !toleranceMet){
    strtRep <- nrow(allBatches)
    # Repeat a batch of simulations
    if (nrow(toRepeat) > 0){
      repeatBatch <- foreach(i=1:nrow(toRepeat)) %dopar% {
        repeatSim(toRepeat[i,], strtRep+i, initializeFunc=initializeScheme, productPipeline=productPipeline, populationImprovement=popImprov1Cyc, targetBudget=targetBudget, bsp=bsp)
      }
    }
    
    strtRep <- strtRep + nrow(toRepeat)
    # Get a new batch of simulations
    newBatch <- foreach(i=1:(batchSize - nrow(toRepeat))) %dopar% {
        runOneRep(strtRep+i, percentRanges=percentRanges, initializeFunc=initializeScheme, productPipeline=productPipeline, populationImprovement=popImprov1Cyc, targetBudget=targetBudget, bsp=bsp)
    }
    
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
    nRepeat <- 0
    toRepeat <- tibble()
    while (nRepeat < nByPareto){
      fitStdErr <- cbind(loPred$fit, loPred$se.fit)
      tst <- returnNonDom(fitStdErr, dir1Low=F, dir2Low=F)
      rows <- which(rownames(fitStdErr) %in% rownames(tst)) 
      rows <- sample(rows, size=min(length(rows), nByPareto - nRepeat))
      fitStdErr <- fitStdErr[-rows,]
      toRepeat <- toRepeat %>% bind_rows(allBatches[as.integer(rownames(tst)),])
    }
    
    allPR <- cbind(allPR, c(unlist(percentRanges), nSimClose=length(bestClose), bestGain=max(loPred$fit), bestSE=bestSE))
    
    # Save batches and results
    if (!is.null(baseFile)){
      saveRDS(allBatches, file=paste0(baseFile, "allBatches.rds"))
      saveRDS(allPR, file=paste0(baseFile, "allPercentRanges.rds"))
    }
    
    batchesDone <- batchesDone + 1
    toleranceMet <- all(percentRanges[,2] - percentRanges[,1] < tolerance)
  }#END keep going until stop instructions
  return(list(allBatches=allBatches, allPercentRanges=allPR))
}
