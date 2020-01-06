#' initFuncSimp function
#'
#' Function to initialize simulation of a breeding program. A single additive trait is simulated. No checks are used in this scheme
#' 
#' @param bsp A list of breeding scheme parameters.  See \code{specifyPipeline} and \code{specifyPopulation} 
#' @return A list containing: 1. The simulation parameters in \code{SP}; 2. The initial records of the breeding program in \code{records}. See \code{fillPipeline} for details; 3. A completed \code{bsp} object
#' 
#' @details Creates the founders and the initial records at the beginning of the simulation of a breeding program.
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#'
#' @export
initFuncSimp <- function(bsp){
  initList <- with(bsp,{
    # Create haplotypes for founder population of outbred individuals
    # Note: default effective population size for runMacs is 100
    founderPop <- runMacs(nInd=nFounders, nChr=nChr, segSites=segSites)
    
    # New global simulation parameters from founder haplotypes
    SP <- SimParam$new(founderPop)
    # Additive and dominance trait architecture
    SP$addTraitA(nQtlPerChr=nQTL, var=genVar)
    # Observed SNPs per chromosome
    SP$addSnpChip(nSNP)
    
    founders <- newPop(founderPop, simParam=SP)
    bsp <- c(bsp, checks=list(NULL))
    
    records <- fillPipeline(founders, bsp, SP)
    
    list(SP=SP, records=records, bsp=bsp)
  })#END with bsp
  return(initList)
}

#' initFuncADChk function
#'
#' function to initialize simulation of a breeding program. A single additive-dominance trait is simulated. Check are used in this scheme
#' 
#' @param bsp A list of breeding scheme parameters.  See \code{specifyPipeline} and \code{specifyPopulation} 
#' @return A list containing: 1. The simulation parameters in \code{SP}; 2. The initial records of the breeding program in \code{records}. See \code{fillPipeline} for details; 3. A completed \code{bsp} object
#' 
#' @details Creates the founders and the initial records at the beginning of the simulation of a breeding program.
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#'
#' @export
initFuncADChk <- function(bsp){
  initList <- with(bsp,{
    # Create haplotypes for founder population of outbred individuals
    # Note: default effective population size for runMacs is 100
    founderPop <- runMacs(nInd=nFounders, nChr=nChr, segSites=segSites)
    
    # New global simulation parameters from founder haplotypes
    SP <- SimParam$new(founderPop)
    # Additive and dominance trait architecture
    SP$addTraitAD(nQtlPerChr=nQTL, var=genVar, meanDD=meanDD, varDD=varDD, useVarA=FALSE)
    # Observed SNPs per chromosome
    SP$addSnpChip(nSNP)
    
    founders <- newPop(founderPop, simParam=SP)
    if (any(nChks > 0)){
      checks <- selectInd(founders, nInd=max(nChks), use="rand", simParam=SP)
    } else checks <- NULL
    bsp <- c(bsp, checks=list(checks))
    
    records <- fillPipeline(founders, bsp, SP)
    
    list(SP=SP, records=records, bsp=bsp)
  })#END with bsp
  return(initList)
}

#' fillPipeline function
#'
#' function to create initial records at the start of a simulation
#'
#' @param founders Pop-class object of the founders of the breeding program
#' @param bsp A list of product pipeline parameters.  It contains nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, stageNames, and nCyclesToKeepRecords.  See \code{runBreedingScheme} for details
#' @return A \code{records} object. A list of lists containing nStages+1 lists. The first list contains the F1 progeny of individuals derived from the \code{parents}. The remaining lists contain one Pop-class object each, with a number of individuals equal to nEntries for that stage. The indivduals have been phenotyped using \code{setPheno}
#' 
#' @details This is a quick and dirty way to create a records object that will be used to simulate breeding schemes
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' records <- with(bsp, {founderPop <- runMacs(nInd=nFounders, nChr=nChr, segSites=segSites)
#' SP <- SimParam$new(founderPop)
#' SP$addTraitA(nQtlPerChr=nQTL, var=genVar)
#' SP$addSnpChip(nSNP)
#' founders <- newPop(founderPop, simParam=SP)
#' bsp <- c(bsp, checks=list(NULL))
#' records <- fillPipeline(founders, bsp, SP)
#' records
#' })
#' 
#' @export
fillPipeline <- function(founders, bsp=NULL, SP){
  records <- with(bsp,{
    
    F1 <- randCross(founders, nCrosses=nCrosses, nProgeny=nProgeny, ignoreGender=T, simParam=SP)
    records <- list(list(F1))
    
    for (stage in 1:nStages){
      entries <- last(records[[stage]])
      
      use <- if_else(stage == 1, "rand", "pheno")
      entries <- selectInd(entries, nInd=nEntries[stage], use=use, simParam=SP)
      # If provided, add checks to the population
      if(!is.null(checks) & nChks[stage] > 0){
        entries <- c(entries, checks[1:nChks[stage]])
      }
      entries <- setPheno(entries, varE=errVars[stage], reps=nReps[stage], simParam=SP)
      
      records <- c(records, list(list(entries)))
    }
    names(records) <- c("F1", stageNames)
    records
  })#END with bsp
  return(records)
}
