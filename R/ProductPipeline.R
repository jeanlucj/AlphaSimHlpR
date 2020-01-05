#' prodPipeSimp function
#'
#' Simple function to advance a simulated breeding product pipeline forward by one generation. No use-input function, no checks. Selection on phenotype
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param selectFunc In this case a dummy function. \code{runBreedingScheme} passes it but this simple function doesn't use it
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has new records created by advancing by a generation
#' 
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
#' 
#' @examples
#' none
#'
#' @export
prodPipeSimp <- function(records, bsp, SP){
  records <- with(bsp,{
    
    for (stage in 1:nStages){
      sourcePop <- last(records[[stage]])
      # Stage 1 coming from untested progeny, so random selection
      use <- if_else(stage == 1, "rand", "pheno")
      entries <- selectInd(sourcePop, nInd=nEntries[stage], use=use, simParam=SP)
      entries <- setPheno(entries, varE=errVars[stage], reps=nReps[stage], simParam=SP)
      records[[stage+1]] <- c(records[[stage+1]], list(entries))
    }
    
    # Remove old records if needed
    if (length(records[[2]]) > nCyclesToKeepRecords) records <- removeOldestCyc(records)
    records
  })#END with bsp
  
  return(records)
}

#' prodPipeFncChk function
#'
#' function to advance a simulated breeding product pipeline forward by one generation. See Gaynor et al. 2017 for the general idea.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of product pipeline parameters.  It contains nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, trialTypeNames, and nCyclesToKeepRecords.  See \code{runBreedingScheme} for details
#' @param selectFunc A function that uses the breeding programs \code{records}, \code{records}, and source population individual ids to generate a vector of selection criterion values that \code{selectInd} will select on. If NULL, \code{selectInd} will select on phenotype. Default is selectFunc=NULL
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has new records created by advancing by a generation
#' 
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
#' 
#' @examples
#' none
#'
#' @export
prodPipeFncChk <- function(records, bsp, SP){
  records <- with(bsp,{
    
    for (stage in 1:nStages){
      sourcePop <- last(records[[stage]])
      # selPipeAdv must be defined in bsp
      entries <- selectInd(sourcePop, nInd=nEntries[stage], trait=selPipeAdv, records=records, ids=sourcePop@id, bsp=bsp, simParam=SP)

      # If provided, add checks to the population
      if(!is.null(checks) & nChks[stage] > 0){
        entries <- c(entries, checks[1:nChks[stage]])
      }
      entries <- setPheno(entries, varE=errVars[stage], reps=nReps[stage], simParam=SP)
      records[[stage+1]] <- c(records[[stage+1]], list(entries))
    }
    
    # Remove old records if needed
    if (length(records[[2]]) > nCyclesToKeepRecords) records <- removeOldestCyc(records)
    records
  })#END with bsp
  
  return(records)
}

#' selectAdvIID function
#'
#' function to select individuals to advance in the product pipeline
#'
#' @param popPheno Matrix of phenotypes of the population being selected. This matrix is provided by \code{selectInd} but need not be used
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param ids Character vector of the ids of the individuals being selected
#' @param bsp A list of product pipeline parameters.  It contains nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, trialTypeNames, and nCyclesToKeepRecords.  See \code{runBreedingScheme} for details
#' @return Real vector of the selection criterion to chose which individuals to advance. \code{selectInd} will sort this vector and choose the individuals
#' @details In deciding which individuals to advance, this function enables you to use all the breeding records rather than just the phenotypes of the population being selected. \code{selectInd} calls this function
#' 
#' @examples
#' none
#'  
#' @export
selectAdvIID <- function(popPheno, records, ids, bsp){
  phenoDF <- framePhenoRec(records, bsp)
  if (!any(ids %in% phenoDF$id)) return(runif(nrow(popPheno)))
  allBLUPs <- iidPhenoEval(phenoDF)
  return(allBLUPs[ids])
}
