#' popImprov1 function
#'
#' Function to improve a simulated breeding population by one cycle. This version takes phenotyped individuals and crosses them to create new F1
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param selectFunc Function to select parents among candidates. Default si the AlfSimHlpR function selectParIID
#' @param SP The AlphaSimR SimParam object
#' @return A records object with a new F1 Pop-class object of progeny coming out of a population improvement scheme
#' 
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' 
#' @examples
#' initList <- initializeFunc()
#' SP <- initList$SP
#' records <- productPipeline(initList$records, bsp=initList$bsp, selectFunc=selectAdvance)
#' records <- populationImprovement(records)
#' recordMeans <- mean_records(records)

#' @export
popImprov1 <- function(records, bsp, selectFunc=selectParIID, SP){
  records <- with(bsp,{
    
    # Make one population to select out of
    bigPop <- mergePops(records[[2]])
    candidates <- bigPop@id
    # Select parents among all individuals
    parents <- bigPop[selectParents(records, candidates, bsp)]
    records[[1]] <- list(randCross(parents, nCrosses=nCrosses, nProgeny=nProgeny, ignoreGender=T, simParam=SP))
    
    records
  })#END with bsp
  
  return(records)
}

#' popImprov2 function
#'
#' Function to improve a simulated breeding population by one cycle. This version does two cycles of predicting F1 individuals and making new F1s
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp List of breeding scheme parameters
#' @param selectFunc Function to select parents among candidates. Default si the AlfSimHlpR function selectParGRM
#' @param SP The AlphaSimR SimParam object
#' @return A records object with a new F1 Pop-class object of progeny coming out of a population improvement scheme
#' 
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' 
#' @examples
#' none

#' @export
popImprov2 <- function(records, bsp, selectFunc=selectParGRM, SP){
  records <- with(bsp,{
    
    for (cycle in 1:2){
      candidates <- records[[1]]@id
      # Select parents among F1
      parents <- records[[1]][selectFunc(records, candidates, bsp)]
      records[[1]] <- list(randCross(parents, nCrosses=nCrosses, nProgeny=nProgeny, ignoreGender=T, simParam=SP))
    }
    
    records
  })#END with bsp
  
  return(records)
}

#' selectParIID function
#'
#' function to select parents among individuals with phenotypes, assuming individual effects are IID
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents
#' @param bsp A list of product pipeline parameters
#' @return Character vector of the ids of the selected individuals
#' @details Accesses all individuals in \code{records} to pick the highest ones
#' @examples none
#' @export
selectParIID <- function(records, candidates, bsp){
  phenoDF <- framePhenoRec(records)
  allBLUPs <- iidPhenoEval(phenoDF, bsp)[candidates]
  nToSelect <- min(bsp$nParents, nrow(allBLUPs))
  ids <- rownames(allBLUPs)[order(allBLUPs, decreasing=T)][1:nToSelect]
  return(ids)
}

#' selectParGRM function
#'
#' function to select parents among individuals with phenotypes, assuming individual effects covary according to a GRM
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents
#' @param bsp A list of product pipeline parameters
#' @return Character vector of the ids of the selected individuals
#' @details Accesses all individuals in \code{records} to pick the highest ones
#' @examples none
#' @export
selectParGRM <- function(records, candidates, bsp){
  phenoDF <- framePhenoRec(records)
  grm <- makeGRM(records)
  allBLUPs <- grmPhenoEval(phenoDF, grm, bsp)[candidates]
  nToSelect <- min(bsp$nParents, nrow(allBLUPs))
  ids <- rownames(allBLUPs)[order(allBLUPs, decreasing=T)][1:nToSelect]
  return(ids)
}
