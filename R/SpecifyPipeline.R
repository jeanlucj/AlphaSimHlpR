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
    nChr <- 2 # Number of chromosomes
    # Population characteristics
    effPopSize <- 100 # Effective size of population generating founders
    quickHaplo <- TRUE # Whether to use AlphaSimR feature to go fast
    segSites <- 20 # Number of segregating sites per chromosome
    nQTL <- 5 # Number of QTL per chromosome
    nSNP <- 5 # Number of observed SNP per chromosome
    genVar <- 40 # Initial genetic variance
    gxyVar <- 15 # Initial genetic x environment variances
    gxlVar <- 10
    gxyxlVar <- 5
    meanDD <- 0.8; varDD <- 0.01 # Mean and variance of dominance degree
    relAA <- 0 # Relative variance that is AxA
    bspNew <- mget(setdiff(ls(), "bspNew"))
    #END no control file
  } else{
    parmNames <- c("nChr", "effPopSize", "quickHaplo", "segSites", "nQTL", "nSNP", "genVar", "gxeVar", "gxyVar", "gxlVar", "gxyxlVar", "meanDD", "varDD", "relAA")
    bspNew <- readControlFile(ctrlFileName, parmNames)
  }
  bsp <- c(bsp, bspNew)
  return(bsp)
}

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
    
    stageToGenotype <- "SDN"
    
    trainingPopCycles <- c(3, 3, 2, 1)
    
    nParents <- 20 # Number of parents in the crossing nursery
    nCrosses <- 20 # Number of crosses entering the pipeline
    nProgeny <- 10 # Number of progeny per cross
    usePolycrossNursery <- FALSE # If TRUE then completely random mating
    nSeeds <- 200 # Number of seeds if usePolycrossNursery
    
    # Don't use optimum contributions in simple default. Define parms in case
    useOptContrib <- FALSE
    nCandOptCont <- 100
    targetEffPopSize <- 30
    # Number of number of entries in each stage
    nEntries <- c(nCrosses*nProgeny, 60, 20, 10)
    nReps <- c(1, 1, 2, 2) # Number of reps used in each stage
    nLocs <- c(1, 2, 2, 3) # Number of locations used in each stage
    # Number of clones that are sent to NCRP at end of product development
    nClonesToNCRP <- 3
    # Number of checks used in each stage
    # Checks are replicated the same as experimental entries
    nChks <- c(2, 1, 1, 1)
    entryToChkRatio <- c(20, 20, 20, 10)
    # Error variances estimated from historical data
    # 200 for SDN is a guess
    errVars <- c(200, 146, 82, 40)
    
    # Use rapid visual selection to move pre-seedlings to SDN
    phenoF1toStage1 <- FALSE
    errVarPreStage1 <- 500
    
    names(nEntries) <- names(nChks) <- names(nReps) <- names(errVars) <- names(trainingPopCycles) <- stageNames
    useCurrentPhenoTrain <- FALSE
    nCyclesToKeepRecords <- 5 # How many cycles to keep records
    nCyclesToRun <- 6 # How many cycles to run the breeding scheme
    # Function to advance individuals from one stage to the next
    selCritPipeAdv <- selCritIID
    selCritPopImprov <- selCritIID
    bspNew <- mget(setdiff(ls(), "bspNew"))
    #END no control file
  } else{
    parmNames <- c("nStages", "stageNames", "stageToGenotype", "trainingPopCycles", "nParents", "nCrosses", "nProgeny", "usePolycrossNursery", "nSeeds", "useOptContrib", "nCandOptCont", "targetEffPopSize", "nEntries", "nReps", "nLocs", "nClonesToNCRP", "nChks", "entryToChkRatio", "errVars", "phenoF1toStage1", "errVarPreStage1", "useCurrentPhenoTrain", "nCyclesToKeepRecords", "nCyclesToRun", "selCritPipeAdv", "selCritPopImprov")
    # Any parameter not specified will have a default set in calcDerivedParms
    bspNew <- readControlFile(ctrlFileName, parmNames)
  }
  bsp <- c(bsp, bspNew)
  bsp <- calcDerivedParms(bsp)
  return(bsp)
}

#' specifyCosts function
#'
#' Function to specify the species and population characteristics. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @param ctrlFileName The name of the text file with parameter values specifying the breeding costs. Must include the path to the file. If NULL a toy example simulation will be set up
#' @return A list containing objects that specify the species and population characteristics.
#'
#' @details Call this function before beginning the simulation
#'
#' @examples
#' bsp <- specifyCosts(bsp)
#'
#' @export
specifyCosts <- function(bsp=NULL, ctrlFileName=NULL){
  if (is.null(ctrlFileName)){ # NULL control file: make toy example
    # Plot costs
    plotCosts <- c(1, 8, 14, 32)
    # Per location cost
    perLocationCost <- 1000
    # Crossing cost
    crossingCost <- 0.2
    # Genotyping cost
    qcGenoCost <- 1.5
    wholeGenomeCost <- 10
    bspNew <- mget(setdiff(ls(), "bspNew"))
    #END no control file
  } else{
    parmNames <- c("plotCosts", "perLocationCost", "crossingCost", "qcGenoCost", "wholeGenomeCost")
    bspNew <- readControlFile(ctrlFileName, parmNames)
  }
  names(bspNew$plotCosts) <- bsp$stageNames
  bsp <- c(bsp, bspNew)
  
  bsp <- calculateBudget(bsp)
  return(bsp)
}

#' calculateBudget function
#'
#' Once the costs are specified in bsp, this function calculates the total annual budget for the breeding scheme
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @return The bsp list with augmented or modified intermediate and total costs
#'
#' @details Call this function once costs have been specified
#'
#' @examples
#' bsp <- calculateBudget(bsp)
#'
#' @export
calculateBudget <- function(bsp){
  # Calculate the program yearly cost
  # Assumptions
  # 1. Genotyping means both QC and whole-genome genotyped, so the cost is
  # nGeno*(qcGenoCost + wholeGenomeCost), where nGeno depends on the stage
  # genotyped.
  # 2. The number of plots in each trial is nEntries*nReps + nChks*chkReps so 
  # the cost of the trial is (nEntries*nReps + nChks*chkReps)*plotCost*nLocs
  # NOTE for develCosts not accounting for number of rapid cycles
  if (length(bsp$plotCosts) != bsp$nStages)
    stop("plotCosts does not have the right length")
  
  bsp$develCosts <- bsp$nCrosses * bsp$nProgeny * bsp$crossingCost
  
  if (is.null(bsp$stageToGenotype) | bsp$stageToGenotype == "F1"){
    nGeno <- bsp$nCrosses * bsp$nProgeny
  } else{
    nGeno <- bsp$nEntries[bsp$stageToGenotype]
  }
  bsp$genotypingCosts <- nGeno * (bsp$qcGenoCost + bsp$wholeGenomeCost)
  
  bsp$trialCosts <- ((bsp$nEntries * bsp$nReps + bsp$nChks * bsp$chkReps) * bsp$nLocs) %*% bsp$plotCost
  
  bsp$locationCosts <- max(bsp$nLocs) * bsp$perLocationCost
  
  bsp$totalCosts <- bsp$develCosts + bsp$genotypingCosts + bsp$trialCosts + bsp$locationCosts
  
  return(bsp)
}

#' calculateChkReps function
#'
#' Once entries are specified, calculate the number of times checks will be replicated. Each rep must have at least one check.
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @return The bsp list with updated chkReps
#'
#' @details Call this function once nEntries, nReps, nChks, and entryToChkRatio have been specified
#'
#' @examples
#' bsp <- calculateChkReps(bsp)
#'
calculateChkReps <- function(bsp){
  # Rerun through this to make sure checks numbers are right
  pairwiseComp <- function(vec1, vec2, fnc){
    return(apply(cbind(vec1, vec2), 1, fnc))
  }
  nChkPlots <- bsp$nEntries * bsp$nReps / bsp$entryToChkRatio
  nChkPlots <- pairwiseComp(nChkPlots, bsp$nReps, max) # At least one check / rep
  chkReps <- ceiling(nChkPlots / bsp$nChks)
  # Safety if nChks or entryToChkRatio misunderstood
  bsp$nChks <- if_else(bsp$entryToChkRatio == 0, 0, bsp$nChks)
  chkReps <- if_else(is.infinite(chkReps) | is.nan(chkReps) | is.na(chkReps), 0, chkReps)
  bsp$chkReps <- chkReps
  
  return(bsp)
}

#' specifyBSP function
#'
#' Function to manually specify a breeding scheme parameters (bsp) object in R, rather than using a control file.
#' Currently does not handle costs. For costs (for now), run specifyCosts() on this.
#' Ideally this is useful for programmatically varying breeding schemes.
#'
#' @param schemeDF data.frame columns: stageNames, nReps, nLocs, nChks, nEntries, entryToChkRatio, errVars
#' @param nParents integer number of parents to cross
#' @param nCrosses integer how many crosses to make
#' @param nProgeny integer how many progeny per cross
#' @param useOptContrib logical whether to use optimal contributions
#' @param nCandOptCont integer how many candidates to consider for opt contrib
#' @param targetEffPopSize numeric target effective population size for OC
#' @param useCurrentPhenoTrain logical whether to use phenotypes for parent sel
#' @param nCyclesToKeepRecords integer eliminate data on cycles above this num
#' @param selCritPipeAdv function used to determine selection criterion for pipe
#' @param selCritPopImprov function used to determine sel crit for pop improv
#' @param nChr integer number of chromosomes for the species
#' @param effPopSize numeric historic effective population size for the species
#' @param segSites integer number of sites segregating per chromosome
#' @param nQTL integer number of loci affecting the trait per chromosome
#' @param nSNP integer number of observed SNPs per chromosome
#' @param genVar numeric genetic variance of the founders
#' @param gxeVar numeric genotype by environment variance of the founders
#' @param meanDD numeric mean dominance deviation. Set to zero for additive
#' @param varDD numeric variance across loci of their dominance deviation
#'
#' @return a named list of of the parameters to specify a breeding scheme simulation
#'
#' @details All arguments are exactly as specified in the control files. Main exception is schemeDF, which is just a tibble() or data.frame version of the set of bsp arguments which are vectors (giving values for each breeding stage). Columns must have names exactly as in the corresponding arguments in control file: stageNames, nReps, nLocs, nChks, nEntries, entryToChkRatio, errVars
#'
#' @examples
#' schemeDF <- tibble(stageNames=c("SDN", "CET", "PYT"),
#'                  nReps=c(1, 1, 2),
#'                  nLocs=c(1, 1, 2),
#'                  nChks=c(1, 1, 2),
#'                  nEntries=c(100, 50, 20),
#'                  entryToChkRatio=c(50, 20, 10),
#'                  errVars=c(150,75,40))
#' bsp <- specifyBSP(schemeDF,nParents = 10, nCrosses = 10, nProgeny = 10,
#'                 useOptContrib = FALSE,
#'                 useCurrentPhenoTrain = TRUE,
#'                 nCyclesToKeepRecords = 1,
#'                 selCritPipeAdv = selCritGRM,
#'                 selCritPopImprov = selCritGRM,
#'                 nChr = 2,effPopSize = 50,
#'                 segSites = 100, nQTL = 5, nSNP = 10,
#'                 genVar = 50, gxeVar = 0, meanDD = 0.05, varDD = 0.25)
#' test <- runBreedingScheme(replication = 1,nCycles = 1,
#'                         initializeFunc = initFuncADChk,
#'                         productPipeline = prodPipeFncChk,
#'                         populationImprovement = popImprov1Cyc,
#'                         bsp = bsp)
#' @export
specifyBSP <- function(schemeDF,
                     nParents,nCrosses,nProgeny,
                     useOptContrib=FALSE,nCandOptCont=NULL,targetEffPopSize=NULL, # if useOptContrib=TRUE, must specify these args
                     useCurrentPhenoTrain=TRUE,
                     nCyclesToKeepRecords,
                     selCritPipeAdv,selCritPopImprov,
                     nChr,effPopSize,
                     segSites,nQTL,nSNP,genVar,gxeVar,meanDD,varDD){
  bspNew <- list()
  bspNew[["nStages"]] <- nrow(schemeDF)
  bspNew[["stageNames"]] <- schemeDF$stageNames
  bspNew[["nReps"]] <- schemeDF$nReps %>% `names <- `(bspNew$stageNames)
  bspNew[["nLocs"]] <- schemeDF$nLocs %>% `names <- `(bspNew$stageNames)
  bspNew[["nChks"]] <- schemeDF$nChks %>% `names <- `(bspNew$stageNames)
  bspNew[["nEntries"]] <- schemeDF$nEntries %>% `names <- `(bspNew$stageNames)
  bspNew[["entryToChkRatio"]] <- schemeDF$entryToChkRatio %>% `names <- `(bspNew$stageNames)
  bspNew[["errVars"]] <- schemeDF$errVars %>% `names <- `(bspNew$stageNames)
  bspNew[["nParents"]] <- nParents
  bspNew[["nCrosses"]] <- nCrosses
  bspNew[["nProgeny"]] <- nProgeny
  bspNew[["useOptContrib"]] <- useOptContrib # if setting this true, there are other arguments that are needed
  bspNew[["useCurrentPhenoTrain"]] <- useCurrentPhenoTrain
  bspNew[["nCyclesToKeepRecords"]] <- nCyclesToKeepRecords
  bspNew[["selCritPipeAdv"]] <- selCritPipeAdv
  bspNew[["selCritPopImprov"]] <- selCritPopImprov
  bspNew[["nChr"]] <- nChr
  bspNew[["effPopSize"]] <- effPopSize
  bspNew[["segSites"]] <- segSites
  bspNew[["nQTL"]] <- nQTL
  bspNew[["nSNP"]] <- nSNP
  bspNew[["genVar"]] <- genVar
  bspNew[["gxeVar"]] <- gxeVar
  bspNew[["meanDD"]] <- meanDD
  bspNew[["varDD"]] <- varDD

  bspNew <- calcDerivedParms(bspNew)
  return(bspNew) 
}

#' calcDerivedParms function
#'
#' Once you have read in parameters from a control file, or set them yourself, there are still a few derived parameters that are needed.  This function calculates them.
#'
#' @param bsp A list. bsp is short for breeding sheme parameters.
#' @return A list bsp that extends the input with a few derived parameters
#'
#' @details This function is only called internally by other functions used to specify the pipeline
#'
#' Should have default if not specified
#' DONE stageToGenotype=stageNames[1]
#' DONE useOptContrib=FALSE, 
#' DONE nCandOptCont=nEntries[1], targetEffPopSize=nParents
#' DONE nChks=0, entryToChkRatio=0
#' DONE phenoF1toStage1=FALSE, errVarPreStage1=genoVar*20
#' DONE quickHaplo=FALSE
#' DONE useCurrentPhenoTrain=FALSE
#' DONE nCyclesToKeepRecords=max(nStages+1, 5)
#' DONE nCyclesToRun=nCyclesToKeepRecords+1
#' DONE selCritPipeAdv=selCritPopImprov=selCritIID
calcDerivedParms <- function(bsp){
  # Function to check if a parameter has no value
  nv <- function(parm){
    is.null(parm) | length(parm) == 0
  }
  
  # Prevent errors having to do with inconsistent parameters
  if (bsp$nSNP + bsp$nQTL >= bsp$segSites){
    print("The number of segregating sites (segSites) has to be greater than the number of SNPs (nSNP) and the number of QTL (nQTL). segSites set 10% bigger than nSNP + nQTL")
    bsp$segSites <- round((bsp$nSNP + bsp$nQTL) * 1.1) + 1
  }
  
  # Some parms have to be logical
  makeLogical <- function(parm){
    if (nv(parm)) parm <- FALSE else parm <- as.logical(parm)
    if (is.na(parm)) parm <- FALSE
    return(parm)
  }
  bsp$useCurrentPhenoTrain <- makeLogical(bsp$useCurrentPhenoTrain)
  bsp$useOptContrib <- makeLogical(bsp$useOptContrib)
  bsp$phenoF1toStage1 <- makeLogical(bsp$phenoF1toStage1)
  bsp$quickHaplo <- makeLogical(bsp$quickHaplo)
  bsp$usePolycrossNursery <- makeLogical(bsp$usePolycrossNursery)
  
  # In case the function is referred by name, replace with actual function
  if (nv(bsp$selCritPipeAdv)) bsp$selCritPipeAdv <- selCritIID
  if (nv(bsp$selCritPopImprov)) bsp$selCritPopImprov <- selCritIID
  if ("character" %in% class(bsp$selCritPipeAdv))
    bsp$selCritPipeAdv <- get(bsp$selCritPipeAdv)
  if ("character" %in% class(bsp$selCritPopImprov))
    bsp$selCritPopImprov <- get(bsp$selCritPopImprov)
  
  # Make sure you keep enough cycles WARNING: not sure this is needed anymore
  bsp$nCyclesToKeepRecords <- max(bsp$nStages+1, bsp$nCyclesToKeepRecords)
  
  if (nv(bsp$nCyclesToRun))
    bsp$nCyclesToRun <- bsp$nCyclesToKeepRecords + 1
  
  # How many clones will go to the National Coordinated Research Program
  if (nv(bsp$nClonesToNCRP)){
    nEndProd <- dplyr::last(bsp$nEntries)
    bsp$nClonesToNCRP <- min(ceiling(nEndProd/2), 3)
  } else{ # Don't specify more clones than there are in the last stage
    bsp$nClonesToNCRP <- min(nEndProd, bsp$nClonesToNCRP)
  }
  
  # If usePolycrossNursery then one seed per cross
  if (nv(bsp$nSeeds)){
    bsp$nSeeds <- bsp$nCrosses * bsp$nProgeny
  }
  if (bsp$usePolycrossNursery){
    bsp$nCrosses <- bsp$nSeeds
    bsp$nProgeny <- 1
  }
  
  # Optimal contributions defaults
  if (bsp$useOptContrib){
    if (bsp$usePolycrossNursery) stop("Polycross nursery and optimal contributions cannot be used together")
    if (nv(bsp$nCandOptCont)) bsp$nCandOptCont <- min(bsp$nEntries[1], bsp$nParents*10)
    if (nv(bsp$targetEffPopSize)) bsp$targetEffPopSize <- bsp$nParents
    # Don't want number of progeny to be too small
    if (bsp$nProgeny < bsp$nSeeds / bsp$targetEffPopSize / 2){
      bsp$nProgeny <- round(bsp$nSeeds / bsp$targetEffPopSize / 2)
      bsp$nCrosses <- round(bsp$nSeeds / bsp$nProgeny)
    }
  }
  
  # Stop and warn user if not enough crosses specified
  if((bsp$nCrosses * bsp$nProgeny) < bsp$nEntries[1]){
    stop("Not enough F1s to fill up Stage 1 trial. [nCrosses * nProgeny >= nEntries for Stage 1] is required")
  }
  
  # Stop and warn user if stageToGenotype is not a named stage
  if (nv(bsp$stageToGenotype)){
    bsp$stageToGenotype <- bsp$stageNames[1]
  }
  if (!(bsp$stageToGenotype %in% c("F1", bsp$stageNames))){
      stop("The stageToGenotype is not one of the pipeline stages")
  }
  
  # Set up trainingPopCycles
  if (nv(bsp$trainingPopCycles)){
    bsp$trainingPopCycles <- integer(bsp$nStages + 1)
    stageNum <- which(bsp$stageNames == bsp$stageToGenotype) + 1
    bsp$trainingPopCycles[stageNum:(bsp$nStages + 1)] <- bsp$nCyclesToKeepRecords
  } else{
    cycF1 <- if_else(bsp$stageToGenotype == "F1", 2, 0)
    bsp$trainingPopCycles <- c(F1=cycF1, bsp$trainingPopCycles)
  }
  
  # Genetic architecture defaults
  if (nv(bsp$meanDD)) bsp$meanDD <- 0
  if (nv(bsp$varDD)) bsp$varDD <- 0
  if (nv(bsp$relAA)) bsp$relAA <- 0
  
  # Figure out how many checks to add to each stage
  pairwiseComp <- function(vec1, vec2, fnc){
    return(apply(cbind(vec1, vec2), 1, fnc))
  }
  if (nv(bsp$entryToChkRatio)) bsp$entryToChkRatio <- integer(bsp$nStages)
  nPlots <- bsp$nEntries * bsp$nReps
  nChkPlots <- nPlots / bsp$entryToChkRatio
  nChkPlots <- pairwiseComp(nChkPlots, bsp$nReps, max) # At least one check / rep
  chkReps <- ceiling(nChkPlots / bsp$nChks)
  # Safety if nChks or entryToChkRatio misunderstood
  bsp$nChks <- if_else(bsp$entryToChkRatio == 0, 0, bsp$nChks)
  chkReps <- if_else(is.infinite(chkReps) | is.nan(chkReps) | is.na(chkReps), 0, chkReps)
  
  # Enforce other defaults
  if (bsp$phenoF1toStage1){
    if (nv(bsp$errVarPreStage1)){
      bsp$errVarPreStage1 <- bsp$errVars[1] * 20
    }
  }
  
  # Check that these vectors are of the right length
  rightLength <- function(vec) length(vec) == bsp$nStages
  v <- list(bsp$stageNames, bsp$nEntries, bsp$entryToChkRatio, bsp$nReps, bsp$nLocs, bsp$errVars)
  names(v) <- c("stageNames", "nEntries", "entryToChkRatio", "nReps", "nLocs", "errVars")
  rl <- sapply(v, rightLength)
  if (any(!rl)){
    stop(paste("These vectors do not have the right length:", paste(names(v)[!rl], collapse=" ")))
  }
  if (length(bsp$trainingPopCycles) != bsp$nStages+1)
    stop("trainingPopCycles does not have the right length")
  
  # Not in use yet...
  if (nv(bsp$analyzeInbreeding)) bsp$analyzeInbreeding <- 0
  
  # Defaults for GxE variance
  if (any(nv(bsp$gxyVar), nv(bsp$gxlVar), nv(bsp$gxyxlVar))){
    if (!nv(bsp$gxeVar)){
      if (nv(bsp$gxyVar)) bsp$gxyVar <- bsp$gxeVar / 3
      if (nv(bsp$gxlVar)) bsp$gxlVar <- bsp$gxeVar / 3
      if (nv(bsp$gxyxlVar)) bsp$gxyxlVar <- bsp$gxeVar / 3
    } else{
      if (nv(bsp$gxyVar)) bsp$gxyVar <- 0
      if (nv(bsp$gxlVar)) bsp$gxlVar <- 0
      if (nv(bsp$gxyxlVar)) bsp$gxyxlVar <- 0
    }
  }

  # Make sure everything has names
  names(bsp$nEntries) <- names(bsp$nChks) <- names(bsp$nReps) <- names(bsp$nLocs) <- names(bsp$errVars) <- names(chkReps) <- names(bsp$entryToChkRatio) <- bsp$stageNames
  names(bsp$trainingPopCycles) <- c("F1", bsp$stageNames)
  bsp <- c(bsp, list(chkReps=chkReps), list(checks=NULL))
  return(bsp)
}

#' adjustEntriesToBudget function
#'
#' Specify a budget, the number of entries for a set of stages, and the stages to adjust to hit the budget. The rules are that the first stage can't be bigger than the number of F1s, and no later stage can be bigger than an earlier stage. The function will not adjust if the rules are broken but will report.
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @param targetBudget Numeric value that you want the budget adjusted to
#' @param fixedEntryStages Named integer vector indicating entry numbers for specific stages. Names must be names of the specific stages
#' @param adjustStages Character vector with names of stages to be changed such that the target budget is achieved
#' 
#' @return A revised bsp with the sizes of the target stages changed to match.
#'
#' @details Call this function after running specifyCosts.
#'
#' @examples
#' bsp <- adjustEntriesToBudget(bsp, targetBudget=50000, fixedEntryStages=c(PYT=100), adjustStages=c("CET", "AYT", "UYT"))
#'
#' @export
adjustEntriesToBudget <- function(bsp, targetBudget, fixedEntryStages=NULL, adjustStages=setdiff(bsp$stageNames, names(fixedEntryStages))){
  if (!is.null(fixedEntryStages)){
    bsp$nEntries[names(fixedEntryStages)] <- fixedEntryStages
  }
  bsp <- calculateChkReps(bsp)
  bsp <- calculateBudget(bsp)
  budgDiff <- (targetBudget - bsp$totalCosts) / length(adjustStages)
  
  for (stage in adjustStages){
    costPerInd <- bsp$nReps[stage] * bsp$nLocs[stage] * bsp$plotCost[stage]
    if (bsp$entryToChkRatio[stage] > 0) costPerInd <- costPerInd * (1 + 1 / bsp$entryToChkRatio[stage])
    if (stage == bsp$stageToGenotype){
        costPerInd <- costPerInd + bsp$qcGenoCost + bsp$wholeGenomeCost
    }
    chngEntries <- round(budgDiff / costPerInd)
    nEntriesNow <- bsp$nEntries[stage] + chngEntries
    if (nEntriesNow < 0) stop(paste("adjustEntriesToBudget: needed nEntries for stage", stage, "below zero"))
    bsp$nEntries[stage] <- nEntriesNow
  }
  
  bsp <- calculateBudget(bsp)
  # Check to make sure no later stages are bigger than earlier stages
  for (stage in 1:bsp$nStages){
    if (stage == 1){
      if (bsp$nEntries[stage] > bsp$nCrosses * bsp$nProgeny) stop("adjustEntriesToBudget: Stage 1 requires more individuals than the number of F1 progeny created")
    } else{
      if (bsp$nEntries[stage] > bsp$nEntries[stage-1]) stop(paste("adjustEntriesToBudget: Stage", stage, " requires more individuals than available from stage", stage - 1))
    }
  }
  return(bsp)
}

#' sampleEntryNumbers function
#'
#' Specify a range of percentages that are allowable for the stages. Function will sample within those percentages and generate a consistent scheme to test. For the stage that gets genotyped, the budget is forced to be the sum of the genotyping + trialling costs.
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @param targetBudget Numeric value that you want the budget adjusted to
#' @param percentRanges Numeric matrix with nStages+1 rows and two columns. Columns are min and max percentage of budget. Rows are costs for crossing, genotyping, trialing each stage.
#' @param nAttempts Integer maximum number of attempts to sample percentages and have them follow the rules of stages becoming progressively smaller.
#' 
#' @return A revised bsp with the sizes of the stages within the percentage ranges specified.
#'
#' @details Call this function after running specifyCosts.
#'
#' @examples
#' Assume stages of CET, PYT, UYT, so percentRanges needs 4 rows
#' Assume CET is genotyped so more budget there
#' percentRanges <- matrix(c(0.02, 0.24, 0.12, 0.12, 0.06, 0.55, 0.30, 0.30), nrow=4)
#' bsp <- sampleEntryNumbers(bsp, targetBudget=50000, percentRanges=percentRanges)
#'
#' @export
sampleEntryNumbers <- function(bsp, targetBudget, percentRanges, nAttempts=5){
  targetBudget <- targetBudget - max(bsp$nLocs) * bsp$perLocationCost
  if (targetBudget < 0) stop("Location costs are above the target budget")
  attemptNo <- 0
  samplingDone <- FALSE
  while (!samplingDone & attemptNo < nAttempts){
    percentages <- apply(percentRanges, 1, function(r) runif(1, r[1], r[2]))
    percentages <- percentages / sum(percentages)
    names(percentages) <- c("F1", bsp$stageNames)
    whchStgGeno <- which(bsp$stageToGenotype == names(percentages)) - 1
    if (length(whchStgGeno) == 0) whchStgGeno <- -1
    
    budgets <- targetBudget * percentages
    
    # How many progeny per cross to make
    f1cost <- bsp$crossingCost
    if (whchStgGeno == 0) f1cost <- f1cost + bsp$qcGenoCost + bsp$wholeGenomeCost
    totProg <- budgets[1] / f1cost
    bsp$nProgeny <- round(totProg / bsp$nCrosses)
    
    for (stage in 1:bsp$nStages){
      costPerInd <- bsp$nReps[stage] * bsp$nLocs[stage] * bsp$plotCost[stage]
      if (bsp$entryToChkRatio[stage] > 0) costPerInd <- costPerInd * (1 + 1 / bsp$entryToChkRatio[stage])
      if (stage == whchStgGeno){
        costPerInd <- costPerInd + bsp$qcGenoCost + bsp$wholeGenomeCost
      }
      bsp$nEntries[stage] <- round(budgets[stage+1] / costPerInd)
    }
    
    # Check to make sure no later stages are bigger than earlier stages
    samplingDone <- TRUE
    for (stage in 1:bsp$nStages){
      if (stage == 1){
        if (bsp$nEntries[stage] > bsp$nCrosses * bsp$nProgeny) samplingDone <- FALSE
      } else{
        if (bsp$nEntries[stage] > bsp$nEntries[stage-1]) samplingDone <- FALSE
      }
    }
    bsp$budgetSamplingDone <- samplingDone
    bsp$budgetPercentages <- percentages
    attemptNo <- attemptNo + 1
  }#END samplingDone
  
  bsp <- calculateBudget(bsp)
  return(bsp)
}

#' function to read a text control file
#'
#' The text file should be organized as follows
#' 1. Any text after a comment symbol # will be ignored
#' 2. Control parameter names should be on their own line
#' 3. Parameter values should be on the following line. If multiple parameter values are needed they should be separated by white space but on the same line
#' @param fileName The name of the text file to be read. Must include the path to the file
#' @param parmNames A string vector with the names of the control parameters that will be searched in the text file
#' @return A named list of the parameter values read from the control file
#'
#' @details Call this function before beginning the simulation
#'
#' @examples
#' params <- readControlFile("./inputDir/ctrlFile.txt", c("nStages", "nParents", "nCrosses"))
#'
#' @export
readControlFile <- function(fileName, parmNames){
  ctrlLines <- readLines(fileName)
  ctrlLines <- sapply(ctrlLines, function(st) strsplit(st, "#", fixed=T)[[1]][1])
  getParm <- function(parmName){
    parmRow <- grep(parmName, ctrlLines)+1
    parms <- unlist(strsplit(ctrlLines[parmRow], "[[:space:]]"))
    names(parms) <- NULL
    parmsNum <- suppressWarnings(as.numeric(parms))
    if (!any(is.na(parmsNum))) parms <- parmsNum
    return(parms)
  }
  parms <- lapply(parmNames, getParm)
  names(parms) <- parmNames
  return(parms)
}
