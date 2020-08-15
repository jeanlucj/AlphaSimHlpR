#' optimizeByLOESS function
#'
#' function to optimize a two-part strategy breeding scheme:
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
#' @param tolerance Numerical difference between min amd max percentage budgets for all stages
#' @param baseDir Directory if you want to have progress saved by batch. Relative  to R working directory. If not empty string, include final /
#' @param maxNumBatches Integer to stop the simulations eventually if the algorithm is not narrowing in on optimal parameter values
#' @param initializeFunc Function for runBreedingScheme
#' @param productPipeline Function for runBreedingScheme
#' @param populationImprovement Function for runBreedingScheme
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
optimizeByLOESS <- function(batchSize, nByPareto=round(batchSize*0.7), targetBudget, percentRanges, startCycle, tolerance, baseDir=NULL, maxNumBatches=10, initializeFunc, productPipeline, populationImprovement, bsp, randomSeed=1234, nCores=1){
  require(parallel)
  
  if (length(randomSeed) == batchSize * maxNumBatches){
    randSeeds <- randomSeed
  } else{
    set.seed(randomSeed)
    randSeeds <- round(runif(batchSize * maxNumBatches, min=-1e9, max=1e9))
    if (!is.null(baseDir)) saveRDS(randSeeds, file=paste0(baseDir, "randSeeds.rds"))
  }
  
  # Guardrails
  nByPareto <- min(batchSize - 1, nByPareto)
  # Sample half by Pareto, and half by probability of beating current best
  nByHiProb <- nByPareto %/% 2
  nByPareto <- nByPareto - nByHiProb

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
      if (!is.null(baseDir)) saveRDS(toRepeat, file=paste0(baseDir, "toRepeat.rds"))
      repeatBatch <- mclapply(1:nrow(toRepeat), function(i) repeatSim(toRepeat[i,], strtRep+i, radius=0.04, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, targetBudget=targetBudget, bsp=bsp, seed=randSeeds[strtRep+i]), mc.cores=nCores)
    } else repeatBatch <- list()
    
    strtRep <- strtRep + nrow(toRepeat)
    # Get a new batch of simulations
    newBatch <- mclapply(1:(batchSize - nrow(toRepeat)), function(i) runOneRep(strtRep+i, percentRanges=percentRanges, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, targetBudget=targetBudget, bsp=bsp, seed=randSeeds[strtRep+i]), mc.cores=nCores)
    
    # Put together with previous simulations
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
    
    # Choose on Pareto frontier of high response and high std. err. of response
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
    # Choose settings that have the best chance of beating the current best
    bestGain <- max(fitStdErr$fit)
    probBetter <- pnorm(bestGain, fitStdErr$fit, fitStdErr$se, lower.tail=FALSE)
    hiProb <- order(probBetter, decreasing=T)[1:nByHiProb]
    toRepeat <- toRepeat %>% bind_rows(allBatches[fitStdErr$batchID[hiProb],])
    
    allPR <- cbind(allPR, c(unlist(percentRanges), nSimClose=length(bestClose), bestGain=max(loPred$fit), bestSE=bestSE))
    
    # Save batches and results
    if (!is.null(baseDir)){
      saveRDS(allBatches, file=paste0(baseDir, "allBatches.rds"))
      saveRDS(allPR, file=paste0(baseDir, "allPercentRanges.rds"))
    }
    
    batchesDone <- batchesDone + 1
    toleranceMet <- all(percentRanges[,2] - percentRanges[,1] < tolerance)
  }#END keep going until maxNumBatches or tolerance
  
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
  on.exit(expr=saveRDS(mget(ls()), file="~/runOneRep.rds"))
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
  rbsOut <- runBreedingScheme(replication=replication, nCycles=bsp$nCyclesToRun, initializeFunc=initializeFunc, productPipeline=productPipeline, populationImprovement=populationImprovement, bsp=bsp)
  
  on.exit()
  return(list(bsp=bsp, stageOutputs=rbsOut$records$stageOutputs))
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
  on.exit(expr=saveRDS(mget(ls()), file="~/repeatSim.rds"))
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
#' 
#' @details Takes a simulation output and returns a vector
#' 
getParmsResponse <- function(oneSim, startCycle){
  on.exit(expr=saveRDS(mget(ls()), file="~/getParmsResp.rds"))
  bsp <- oneSim$bsp
  parms <- c(budget=bsp$budgetPercentages, nProgeny=bsp$nProgeny, bsp$nEntries, totCost=bsp$totalCosts)
  so <- oneSim$stageOutputs
  resp <- (dplyr::filter(so, stage=="F1" & year==bsp$nCyclesToRun) %>% dplyr::select(genValMean, genValSD)) - (dplyr::filter(so, stage=="F1" & year==startCycle) %>% dplyr::select(genValMean, genValSD))
  on.exit()
  return(unlist(c(parms, resp)))
}

