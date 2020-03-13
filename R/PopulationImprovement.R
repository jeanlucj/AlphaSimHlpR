#' popImprov1Cyc function
#'
#' Function to improve a simulated breeding population by one cycle. This version takes phenotyped individuals and crosses them to create new F1
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP The AlphaSimR SimParam object
#' @return A records object with a new F1 Pop-class object of progeny coming out of a population improvement scheme
#' 
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
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
popImprov1Cyc <- function(records, bsp, SP){
  candidates <- records[[1]]@id
  # Select parents among all individuals
  parents <- records[[1]][selectParIID(records, candidates, bsp)]
  progeny <- randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreGender=T, simParam=SP)
  records[[1]] <- c(records[[1]], progeny)
  return(records)
}

#' popImprov2Cyc function
#'
#' Function to improve a simulated breeding population by one cycle. This version does two cycles of predicting F1 individuals and making new F1s
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp List of breeding scheme parameters
#' @param SP The AlphaSimR SimParam object
#' @return A records object with a new F1 Pop-class object of progeny coming out of a population improvement scheme
#' 
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- prodPipeSimp(records, bsp, SP)
#' records <- popImprov2(records, bsp, SP)
#' 
#' @export
popImprov2Cyc <- function(records, bsp, SP){
  for (cycle in 1:2){
    candidates <- records[[1]]@id
    parents <- records[[1]][selectParGRM(records, candidates, bsp, SP)]
    progeny <- randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreGender=T, simParam=SP)
    records[[1]] <- c(records[[1]], progeny)
  }
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
#' @details Accesses all individuals in \code{records} to pick the highest among candidates. If candidates do not have records, a random sample is returned
#' 
#' @examples
#' allPop <- mergePops(records[[2]])
#' candidates <- allPop@id
#' parents <- allPop[selectParIID(records, candidates, bsp)]
#' 
#' @export
selectParIID <- function(records, candidates, bsp){
  phenoDF <- framePhenoRec(records)
  if (!any(candidates %in% phenoDF$id)) return(sample(candidates, bsp$nParents))
  allBLUPs <- iidPhenoEval(phenoDF)
  allBLUPs <- allBLUPs[names(allBLUPs) %in% candidates]
  nToSelect <- min(bsp$nParents, length(allBLUPs))
  ids <- names(allBLUPs)[order(allBLUPs, decreasing=T)][1:nToSelect]
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
#' @examples 
#' candidates <- records[[1]][[1]]@id
#' parents <- records[[1]][[1]][selectParGRM(records, candidates, bsp, SP)]
#' 
#' @export
selectParGRM <- function(records, candidates, bsp, SP){
  phenoDF <- framePhenoRec(records)
  grm <- makeGRM(records, SP)
  if (!any(candidates %in% rownames(grm))) return(sample(candidates, bsp$nParents))
  # Remove individuals with phenotypes but who no longer have geno records
  # I am not sure this can happen, but it is a safeguard
  phenoDF <- phenoDF[phenoDF$id %in% rownames(grm),]
  allBLUPs <- grmPhenoEval(phenoDF, grm)
  allBLUPs <- allBLUPs[names(allBLUPs) %in% candidates]
  nToSelect <- min(bsp$nParents, length(allBLUPs))
  ids <- names(allBLUPs)[order(allBLUPs, decreasing=T)][1:nToSelect]
  return(ids)
}
