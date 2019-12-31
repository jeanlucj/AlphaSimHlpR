#' fillPipeline function
#'
#' function to (do something)
#'
#' @param parents [value]
#' @param checks [value]. Default is checks=NULL
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
fillPipeline <- function(parents, checks=NULL){
  records <- list(list(randCross(parents, nCrosses=nCrosses, nProgeny=nProgeny, ignoreGender=T)))
  
  for (stage in 1:nStages){
    entries <- last(records[[stage]])
    entries <- setPheno(entries, varE=errVars[stage], reps=nReps[stage])
    entries <- selectInd(entries, nInd=nEntries[stage])
    # If provided, add checks to the population
    if(!is.null(checks) & nChks[stage] > 0){
      entries <- c(entries, checks[1:nChks[stage]])
    }
    records <- c(records, list(list(entries)))
  }
  names(records) <- c("F1", trialTypeNames)
  
  return(records)
}
#' productPipeline function
#'
#' function to (do something)
#'
#' @param records [value]
#' @param checks [value]. Default is checks=NULL
#' @param selectFunc [value]. Default is selectFunc=NULL
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
productPipeline <- function(records, checks=NULL, selectFunc=NULL){
  
  for (stage in 1:nStages){
    sourcePop <- last(records[[stage]])
    # selectInd can't deal with missing phenotypes
    if (any(is.na(pheno(sourcePop)))){
      rownames(sourcePop@pheno) <- sourcePop@id
      sourcePop@pheno[is.na(sourcePop@pheno)] <- 0
    }
    
    # Did user provide function to select individuals to advance?
    if (is.null(selectFunc)){ # Just phenotypic selection
      # If there are checks, remove before selecting
      # Only relevant if sourcePop from SDN or later (stage > 1)
      if (stage > 1){
        if(!is.null(checks) & nChks[stage-1] > 0){
          sourcePop <- sourcePop[1:nEntries[stage-1]]
        }
      }
      entries <- selectInd(sourcePop, nInd=nEntries[stage])
    } else{ 
      # User provided a function to select individuals to advance
      # If there are checks, selectFunc will deal with them
      entries <- selectInd(sourcePop, nInd=nEntries[stage], trait=selectFunc, records=records, checks=checks, ids=sourcePop@id)
    }
    # If provided, add checks to the population
    if(!is.null(checks) & nChks[stage] > 0){
      entries <- c(entries, checks[1:nChks[stage]])
    }
    entries <- setPheno(entries, varE=errVars[stage], reps=nReps[stage])
    records[[stage+1]] <- c(records[[stage+1]], list(entries))
  }
  return(records)
}
#' selectAdvance function
#'
#' function to (do something)
#'
#' @param popPheno [value]
#' @param records [value]
#' @param checks [value]
#' @param ids [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
selectAdvance <- function(popPheno, records, checks, ids){
  phenoDF <- framePhenoRec(records)
  if (sum(ids %in% phenoDF$id) == 0) return(runif(nrow(popPheno)))
  # phenoDF <- phenoDF[phenoDF$id %in% ids,]
  allBLUPs <- iidPhenoEval(phenoDF)
  return(allBLUPs[ids, 1])
}
