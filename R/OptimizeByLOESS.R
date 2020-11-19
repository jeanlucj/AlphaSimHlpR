#' optimizeByLOESS function
#'
#' Function to optimize a two-part strategy breeding scheme:
#' 1. Simulate a batch using given percentage ranges
#' 2. Perform LOESS fit to the gains
#' 3. Find budget with best estimated gain
#' 4. Calculate new percentage ranges: any simulation within 2*StdErr of best
#' 5. Decide on some simulations to repeat:
#'    1. Parameter space with high gain and high std err: need more info there
#'    2. Parameter space with high gain: high probability that it's best
#' Go back to 1.
#'
#' @param batchSize Integer number of simulations between LOESS fits
#' @param targetBudget Numeric value that you want the budget adjusted to
#' @param percentRanges Numeric matrix (nStages + 1 rows and 2 columns) of percentage budget allocation to crossing (F1) and all of the stages.  If there is a stage that is genotyped, the genotyping cost is added to that stage
#' @param startCycle Integer the start cycle from which to measure gain. The end cycle will be the last cycle. Set the startCycle so there is enough burn-in
#' @param wgtPopImprov Two components contribute to the utility: yearly population improvement and yearly gain in the product development. This is the weight for the first. The other has weight 1 - wgtPopImprov
#' @param tolerance Numerical difference between min amd max percentage budgets for all stages
#' @param baseDir Directory if you want to have progress saved by batch. Relative  to R working directory. If not empty string, include final /
#' @param maxNumBatches Integer to stop the simulations eventually if the algorithm is not narrowing in on optimal parameter values
#' @param initializeFunc Function for runBreedingScheme
#' @param productPipeline Function for runBreedingScheme
#' @param populationImprovement Function for runBreedingScheme
#' @param bsp  A list of breeding scheme parameters.
#' @param randomSeed Integer seed for random number generator
#' @param nCores Integer number of cores to use for parallel simulation
#' @return Numeric matix with all simulations budget allocations, gen mean change, gen std dev change, total cost.
#' 
#' @details A wrapper to repeatedly simulate a scheme with different budget allocations to find optimal allocations
#' 
#' @examples
#' 
#' @export
optimizeByLOESS <- function(batchSize, targetBudget, percentRanges, startCycle, wgtPopImprov, tolerance, baseDir=NULL, maxNumBatches=10, initializeFunc, productPipeline, populationImprovement, bsp, randomSeed=NULL, nCores=1){
  on.exit(expr={print(traceback()); saveRDS(mget(ls()), file="~/optimizeByLOESS.rds")})
  require(parallel)
  
  if (length(randomSeed) == batchSize * maxNumBatches){
    randSeeds <- randomSeed
  } else{
    if (is.null(randomSeed)) randomSeed <- round(runif(1)*1e9)
    set.seed(randomSeed)
    randSeeds <- round(runif(batchSize * maxNumBatches, min=-1e9, max=1e9))
    if (!is.null(baseDir)) saveRDS(randSeeds, file=paste0(baseDir, "randSeeds.rds"))
  }
  
  # Guardrails
  if (wgtPopImprov < 0 | wgtPopImprov > 1) wgtPopImprov <- 0.5
  
  # Implement cockamamie scheme to sample different numbers of repeats
  nuts <- tibble(ratio=6:16, rand=0.6-0.05*0:10, hiGain=0.085*0:10, pareto=0.2-0.015*0:10, hiStEr=0.2-0.02*0:10)
  
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
      repeatBatch <- mclapply(1:nrow(toRepeat), function(i) repeatSim(toRepeat[i,], strtRep+i, radius=0.04, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, targetBudget=targetBudget, bsp=bsp, seed=randSeeds[strtRep+i]), mc.cores=nCores)
    } else repeatBatch <- list()
    
    strtRep <- strtRep + nrow(toRepeat)
    # Get a new batch of simulations
    newBatch <- mclapply(1:(batchSize - nrow(toRepeat)), function(i) runOneRep(strtRep+i, percentRanges=percentRanges, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, targetBudget=targetBudget, bsp=bsp, seed=randSeeds[strtRep+i]), mc.cores=nCores)
    
    # Put together with previous simulations
    newBatch <- as_tibble(t(sapply(c(repeatBatch, newBatch), getParmsResponse, startCycle=startCycle, wgtPopImprov=wgtPopImprov)), .name_repair="universal")
    allBatches <- allBatches %>% bind_rows(newBatch)
    
    # Non-Parametric LOESS response
    loFormula <- paste0("response ~ ", paste0(bsp$stageNames, collapse=" + "))
    loFM <- loess(loFormula, data=allBatches)
    loPred <- predict(loFM, se=T)
    whichBest <- which.max(loPred$fit)
    bestSE <- loPred$se.fit[whichBest]
    bestClose <- which(max(loPred$fit) - loPred$fit < 2*bestSE)
    percentRanges <- t(apply(allBatches[bestClose, ] %>% dplyr::select(contains("budget")), 2, range))
    
    # Cockamamie
    ratio <- round(diff(range(loPred$fit)) / mean(loPred$se.fit))
    percRow <- which(nuts$ratio == ratio)
    if (ratio < min(nuts$ratio)) percRow <- 1
    if (ratio > max(nuts$ratio)) percRow <- 11
    nByPareto <- round(nuts$pareto[percRow] * batchSize)
    nByHiGain <- round(nuts$hiGain[percRow] * batchSize)
    nByStdErr <- round(nuts$hiStEr[percRow] * batchSize)
    
    # Choose settings to repeat to get more information there
    toRepeat <- tibble()
    fitStdErr <- tibble(batchID=1:nrow(allBatches), fit=loPred$fit, se=loPred$se.fit)
    # Repeat settings on Pareto frontier
    rows <- returnNonDom(fitStdErr, dir1Low=F, dir2Low=F, var1name="fit", var2name="se")$batchID
    # If too many rows, keep only those with the highest gain
    if (length(rows) > nByPareto){
      rows <- rows[order(allBatches$response[rows], decreasing=T)[1:nByPareto]]
    }
    toRepeat <- toRepeat %>% bind_rows(allBatches[rows,])
    nByPareto <- nByPareto - length(rows)
    
    # Repeat settings that have the best chance of beating the current best
    bestGain <- max(fitStdErr$fit)
    probBetter <- pnorm(bestGain, fitStdErr$fit, fitStdErr$se, lower.tail=FALSE)
    hiProb <- order(probBetter, decreasing=T)[1:nByHiGain]
    toRepeat <- toRepeat %>% bind_rows(allBatches[fitStdErr$batchID[hiProb],])
    fitStdErr <- fitStdErr %>% dplyr::filter(!(batchID %in% rows))
    
    # Some more with high Pareto ranks
    while (nByPareto > 0){
      rows <- returnNonDom(fitStdErr, dir1Low=F, dir2Low=F, var1name="fit", var2name="se")$batchID
      if (length(rows) > nByPareto){
        rows <- rows[order(allBatches$response[rows], decreasing=T)[1:nByPareto]]
      }
      toRepeat <- toRepeat %>% bind_rows(allBatches[rows,])
      fitStdErr <- fitStdErr %>% dplyr::filter(!(batchID %in% rows))
      nByPareto <- nByPareto - length(rows)
    }
    
    # Some with high Std Error: we need more information there
    if (nByStdErr > 0){
      rows <- fitStdErr$batchID[order(fitStdErr$se, decreasing=T)[1:nByStdErr]]
    }
    toRepeat <- toRepeat %>% bind_rows(allBatches[rows,])

    allPR <- cbind(allPR, c(unlist(percentRanges), nSimClose=length(bestClose), bestGain=max(loPred$fit), bestSE=bestSE))
    
    batchesDone <- batchesDone + 1

        # Save batches and results
    if (!is.null(baseDir)){
      saveRDS(allBatches, file=paste0(baseDir, "allBatches.rds"))
      saveRDS(allPR, file=paste0(baseDir, "allPercentRanges.rds"))
      saveRDS(toRepeat, file=paste0(baseDir, paste0("toRepeat", batchesDone, ".rds")))
    }
    
    toleranceMet <- all(percentRanges[,2] - percentRanges[,1] < tolerance)
  }#END keep going until maxNumBatches or tolerance
  
  on.exit()
  return(list(allBatches=allBatches, allPercentRanges=allPR))
}

#' runOneRep function
#'
#' Utility function for optimizeByLOESS
#' Sample a budget according to allowable ranges for the different stages and run a rep of the breeding scheme with that budget
#'
#' @param replication Integer replication
#' @param percentRanges Real matrix with F1 + nStages allowable budget ranges
#' @param initializeFunc Function for runBreedingScheme
#' @param productPipeline Function for runBreedingScheme
#' @param populationImprovement Function for runBreedingScheme
#' @param targetBudget Real how much to spend total. See CostControlFile
#' @param bsp List breeding scheme parameters
#' @param seed=NULL Integer seed to initialize pseudo-random number generator
#' 
#' @details Allows sampling the budget and running the simulation in one function
#' 
runOneRep <- function(replication, percentRanges, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp, seed=NULL){
  on.exit(expr={print(traceback()); saveRDS(mget(ls()), file="~/runOneRep.rds")})
  if (!is.null(seed)) set.seed(seed)
  bsp$budgetSamplingDone <- FALSE
  while (!bsp$budgetSamplingDone){
    bsp <- sampleEntryNumbers(bsp, targetBudget, percentRanges)
    if (!bsp$budgetSamplingDone){
      # Sampling failed, so shift some budget from later stages to earlier stages
      for (i in nrow(percentRanges):2){
        for (j in 1:2){
          if (percentRanges[i, j] > percentRanges[i-1, j]){
            percentRanges[i, j] <- max(0.01, percentRanges[i, j] - 0.01)
            percentRanges[i-1, j] <- min(0.99, percentRanges[i-1, j] + 0.01)
          }
        }
      }
    }#END make sure proper sampling of budgets was done
  }# Carry on
  # Dealing with errors
  options(try.outFile = stdout())
  nTry <- 0
  badRBS <- TRUE
  while (badRBS){
    rbsOut <- try(runBreedingScheme(replication=replication, nCycles=bsp$nCyclesToRun, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, bsp=bsp))
    badRBS <- "try-error" %in% class(rbsOut)
    nTry <- nTry + 1
    cat("\n nTry", nTry, "\n")
    if (nTry > 10) stop("Too many runBreedingScheme tries")
  }
  
  on.exit()
  return(list(bsp=bsp, stageOutputs=rbsOut$records$stageOutputs, nTry=nTry))
}#END runOneRep

#' repeatSim function
#'
#' Utility function for optimizeByLOESS
#' Run a simulation in the vicinity of a previous simulation
#'
#' @param parmRow Vector resulting from getParmsResponse and that contains parameters for the simulation to be repeated
#' @param replication Integer replication
#' @param radius Real radius around which to sample the budget so it's not an _exact_ repeat
#' @param initializeFunc Function for runBreedingScheme
#' @param productPipeline Function for runBreedingScheme
#' @param populationImprovement Function for runBreedingScheme
#' @param targetBudget Real how much to spend total. See CostControlFile
#' @param bsp List breeding scheme parameters
#' @param seed=NULL Integer seed to initialize pseudo-random number generator
#' 
#' @details Sample the budget and run a simulation close to a previous one. Useful to rerun where the prediction is particularly high (we need to know if that prediction is for real) or where the standard error is particularly high (we need more info in that vicinity)
#' 
repeatSim <- function(parmRow, replication, radius=0.02, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp, seed=NULL){
  on.exit(expr={print(traceback()); saveRDS(mget(ls()), file="~/repeatSim.rds")})
  budg <- parmRow %>% dplyr::select(contains("budget"))
  percentRanges <- t(sapply(unlist(budg), function(prc) c(max(0, prc - radius), min(1, prc + radius))))
  rorOut <- runOneRep(replication, percentRanges, initializeFunc, productPipeline, populationImprovement, targetBudget, bsp, seed)
  on.exit()
  return(rorOut)
}

#' getParmsResponse function
#'
#' Utility function for optimizeByLOESS
#' Pull interesting parameters from from the stageOutputs that will serve to calculate optimality of the simulation. Currently the mean and genetic std dev of the F1 population
#'
#' @param oneSim Output of runBreedingScheme
#' @param startCycle Integer first cycle after burn-in is done
#' @param wgtPopImprov Two components contribute to the utility: yearly population improvement and yearly gain in the product development. This is the weight for the first. The other has weight 1 - wgtPopImprov
#' 
#' @details Takes a simulation output and returns a vector
#' 
getParmsResponse <- function(oneSim, startCycle, wgtPopImprov=1){
  on.exit(expr={print(traceback()); saveRDS(mget(ls()), file="~/getParmsResp.rds")})
  bsp <- oneSim$bsp
  parms <- c(budget=bsp$budgetPercentages, nProgeny=bsp$nProgeny, bsp$nEntries, totCost=bsp$totalCosts)
  so <- oneSim$stageOutputs
  
  # Gain across versus within cycles
  soLast <- so %>% dplyr::filter(stage==last(bsp$stageNames) & cycle >= startCycle)
  soF1 <- so %>% dplyr::filter(stage=="F1" & cycle > startCycle & cycle <= max(soLast$cycle)+1)
  popImpr <- mean(dplyr::lead(soF1$genValMean) - soF1$genValMean, na.rm=T)
  prodDev <- mean((soLast$gvOfBestCrit - soF1$genValMean) / bsp$nStage, na.rm=T)
  resp <- wgtPopImprov * popImpr + (1 - wgtPopImprov) * prodDev
  
  on.exit()
  return(unlist(c(parms, response=resp)))
}
