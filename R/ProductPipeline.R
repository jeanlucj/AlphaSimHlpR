#' productPipeline function
#'
#' function to advance a simulated breeding product pipeline forward by one generation. See Gaynor et al. 2017 for the general idea.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of product pipeline parameters.  It contains nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, trialTypeNames, and nCyclesToKeepRecords.  See \code{runBreedingScheme} for details
#' @param selectFunc A function that uses the breeding programs \code{records}, \code{records}, and source population individual ids to generate a vector of selection criterion values that \code{selectInd} will select on. If NULL, \code{selectInd} will select on phenotype. Default is selectFunc=NULL
#' @return A records object that has new records created by advancing by a generation
#' 
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
#' 
#' @examples
#' initList <- initializeFunc()
#' SP <- initList$SP
#' records <- productPipeline(initList$records, bsp=initList$bsp, selectFunc=selectAdvance)
#' records <- populationImprovement(records)
#' recordMeans <- mean_records(records)
#'
#' @export
productPipeline <- function(records, bsp=NULL, selectFunc=NULL, SP){
  records <- with(bsp,{

    for (stage in 1:nStages){
      sourcePop <- last(records[[stage]])
      
      # Did user provide function to select individuals to advance?
      if (is.null(selectFunc)){ # No function so just phenotypic selection
        if (stage == 1){ # The initial stage coming from untested progeny
          use <- "rand"
        } else{
          use <- "pheno"
          # Checks may be among records. If so remove before selecting
          if(!is.null(checks) & nChks[stage-1] > 0){
            sourcePop <- sourcePop[1:nEntries[stage-1]]
          }
        }
        entries <- selectInd(sourcePop, nInd=nEntries[stage], use=use, simParam=SP)
      } else{ 
        # User provided a function to select individuals to advance
        # If there are checks, selectFunc will deal with them
        entries <- selectInd(sourcePop, nInd=nEntries[stage], trait=selectFunc, records=records, ids=sourcePop@id, bsp=bsp, simParam=SP)
      }
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

#' selectAdvance function
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
#' founderPop <- quickHaplo(nInd=50, nChr=1, segSites=100)
#' SP <- SimParam$new(founderPop)
#' SP$addTraitA(nQtlPerChr=10)
#' parents <- newPop(founderPop)
#' records <- fillPipeline(parents)
#' sourcePop <- last(records[[2]])
#' toAdvance <- selectInd(sourcePop, nInd=ceiling(nInd(sourcePop)/2), trait=selectAdvance, records=records, ids=sourcePop@id)
#' 
#' @export
selectAdvance <- function(popPheno, records, ids, bsp){
  phenoDF <- framePhenoRec(records)
  if (!any(ids %in% phenoDF$id)) return(runif(nrow(popPheno)))
  allBLUPs <- iidPhenoEval(phenoDF, bsp)
  return(allBLUPs[ids, 1])
}
