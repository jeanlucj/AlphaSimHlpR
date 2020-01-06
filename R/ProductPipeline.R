#' prodPipeSimp function
#'
#' Simple function to advance a simulated breeding product pipeline forward by one generation. No use-input function, no checks. Selection on phenotype
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has new records created by advancing by a generation
#' 
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- prodPipeSimp(records, bsp, SP)
#' records <- popImprov1(records, bsp, SP)
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
#' @param bsp A list of breeding scheme parameters
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has new records created by advancing by a generation
#' 
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- prodPipeFncChk(records, bsp, SP)
#' records <- popImprov1(records, bsp, SP)
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
#' @param bsp A list of breeding scheme parameters.
#' @return Real vector of the selection criterion to chose which individuals to advance. \code{selectInd} will sort this vector and choose the individuals
#' @details In deciding which individuals to advance, this function enables you to use all the breeding records rather than just the phenotypes of the population being selected. \code{selectInd} calls this function
#' 
#' @examples
#' stage <- 2
#' entries <- with(bsp, {
#' sourcePop <- last(records[[stage]])
#' entries <- selectInd(sourcePop, nInd=nEntries[stage], trait=selectAdvIID, records=records, ids=sourcePop@id, bsp=bsp, simParam=SP)
#' entries
#' })
#'  
#' @export
selectAdvIID <- function(popPheno, records, ids, bsp){
  phenoDF <- framePhenoRec(records, bsp)
  if (!any(ids %in% phenoDF$id)) return(runif(nrow(popPheno)))
  allBLUPs <- iidPhenoEval(phenoDF)
  return(allBLUPs[ids])
}
