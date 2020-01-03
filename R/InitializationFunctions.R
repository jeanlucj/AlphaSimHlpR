#' initializeFunc function
#'
#' function to initialize simulation of a breeding program. A single additive-dominance trait is simulated.
#' 
#' @param spp A list of species parameters. It contains nFounders, nChr, segSites, nQTL, genVar, meanDD, varDD, nSNP.  See \code{runBreedingScheme} for details
#' @param ppp A list of product pipeline parameters.  It contains nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, trialTypeNames, and nCyclesToKeepRecords.  See \code{runBreedingScheme} for details
#' @return A list containing: 1. The simulation parameters in \code{SP}; 2. The initial records of the breeding program in \code{records}. See \code{fillPipeline} for details; 3. A list of product pipeline parameters containing nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, trialTypeNames, and nCyclesToKeepRecords.
#' 
#' @details Creates the founders and the initial records at the beginning of the simulation of a breeding program.
#' 
#' @examples
#' initList <- initializeFunc()
#' SP <- initList$SP
#' records <- productPipeline(initList$records, ppp=initList$ppp)
#' records <- populationImprovement(initList$records, ppp=initList$ppp)
#' recordMeans <- mean_records(records)
#'
#' @export
initializeFunc <- function(spp=NULL, ppp=NULL){
  if (is.null(spp)) spp <- list(nFounders=50, nChr=10, segSites=200, nQTL=50, genVar=1, meanDD=0.2, varDD=0.1, nSNP=50)
  if (is.null(ppp)) ppp <- list(nParents=30, nCrosses=20, nProgeny=10, nStages=2, errVars=c(4, 2), nReps=c(1, 2), nEntries=nCrosses*nProgeny*c(1, 0.5), nChks=c(0, 0), trialTypeNames=c("PYT", "AYT"), nCyclesToKeepRecords=7)
  
  initList <- with(spp,{
    with(ppp,{
      # Create haplotypes for founder population of outbred individuals
      # Default effective population size for runMacs is 100
      founderPop <- runMacs(nInd=nFounders, nChr=nChr, segSites=segSites)
      
      # New global simulation parameters from founder haplotypes
      SP <- SimParam$new(founderPop)
      # Additive and dominance trait architecture
      SP$addTraitAD(nQtlPerChr=nQTL, var=genVar, meanDD=meanDD, varDD=varDD, useVarA=FALSE)
      # Observed SNPs per chromosome
      SP$addSnpChip(nSNP)
      
      parents <- newPop(founderPop, simParam=SP)
      if (any(nChks > 0)){
        checks <- selectInd(parents, nInd=max(nChks), use="rand", simParam=SP)
      } else checks <- NULL
      ppp <- c(ppp, checks=list(checks))
      
      records <- fillPipeline(parents, ppp, SP)
      
      list(SP=SP, records=records, ppp=ppp)
    })#END with ppp
  })#END with spp
  return(initList)
}

#' fillPipeline function
#'
#' function to create initial records at the start of a simulation
#'
#' @param parents Pop-class object of the founders of the breeding program
#' @param ppp A list of product pipeline parameters.  It contains nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, trialTypeNames, and nCyclesToKeepRecords.  See \code{runBreedingScheme} for details
#' @return A \code{records} object. A list of lists containing nStages+1 lists. The first list contains the F1 progeny of individuals derived from the \code{parents}. The remaining lists contain one Pop-class object each, with a number of individuals equal to nEntries for that stage. The indivduals have been phenotyped using \code{setPheno}
#' 
#' @details This is a quick and dirty way to create a records object that will be used to simulate breeding schemes
#' 
#' @examples
#' founderPop <- quickHaplo(nInd=50, nChr=1, segSites=100)
#' SP <- SimParam$new(founderPop)
#' SP$addTraitA(nQtlPerChr=10)
#' parents <- newPop(founderPop)
#' records <- fillPipeline(parents)
#'
#' @export
fillPipeline <- function(parents, ppp=NULL, SP){
  if (is.null(ppp)) ppp <- list(nParents=30, nCrosses=20, nProgeny=10, nStages=2, errVars=c(4, 2), nReps=c(1, 2), nEntries=nCrosses*nProgeny*c(1, 0.5), nChks=c(0, 0), trialTypeNames=c("PYT", "AYT"), nCyclesToKeepRecords=7, checks=NULL)
  
  records <- with(ppp,{

    F1 <- randCross(parents, nCrosses=nCrosses, nProgeny=nProgeny, ignoreGender=T, simParam=SP)
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
    names(records) <- c("F1", trialTypeNames)
    records
  })#END with ppp
  return(records)
}
